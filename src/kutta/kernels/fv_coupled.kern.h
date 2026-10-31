#include "kutta/kernels/prelude.h"

// Pressure-based coupled solve (Fluent Theory Guide 23.4.3.2, eqs. 23-58 to
// 23-64): one block system in cell unknowns (u, v, w, p~) per outer iteration.
// Block vectors are one buffer of 4n floats: [u | v | w | p~], p~ = p/wc.
//
// Delta form: the RHS is the steady SIMPLEC residual itself, so the converged
// fixed point is SIMPLEC's discretization whatever the Jacobian approximates.
// Jacobian, per face f between cells i and j, s = +1 out of i:
//   momentum  a_ij = -mu_f D_f + min(s F_f, 0)       (the assembled transport row)
//             G_ij = wc (0.5 s S_f + Q_i d_ij)       (implicit pressure force)
//   volume continuity H_ij = wc s S_f a_j/(a_i + a_j), C_ij = -wc^2 rc_f D_f
// H follows the residual's momentum-weighted face velocity (Theory Guide,
// Discretization of the Continuity Equation; eq. 18.4-6 in Fluent 12). Equal
// halves misstate a cell's own flux where a_j/(a_i + a_j) is far from 1/2:
// separated and wake cells next to fast ones. The momentum diagonal is
// floored and relaxed by the Courant number (coupled_copy_diagonal).
// The residual's face pressures extrapolate with the cell's least-squares
// gradient, g_i = M_i^-1 sum_j d_ij (p_j - p_i). Q_i = T_i M_i^-1 carries that
// own-cell part exactly, with T_i = sum_f S_f (x_f - x_i)^T over boundary faces
// with extrapolated pressure and half of it over internal faces. Without it,
// tetrahedral boundary cells respond up to ~3x more strongly than the Jacobian
// predicts and the outer iteration diverges. The neighbour-side gradients are
// outside the compact stencil and remain explicit.
// The continuity rows are scaled by wc = rho*U and p by 1/wc, so both rows
// carry momentum units.
//
// Smoothers solve each cell's arrow block with the pressure diagonal raised
// by a bound on the neighbour part of the Schur complement C - H D^-1 G:
//   cs_i = sum_{k != i} sum_c |H_ik,c| R_k,c / D_k,c,  R_k,c = sum_j |G_kj,c|.
// The u-p coupling is near-skew (G_ji = -H_ij for equal face weights), so
// damped block Jacobi converges only while
// rho(c^-1 H D^-1 G) < (2 - omega)/omega. With the compact Rhie-Chow diagonal
// alone that fails on stretched cells at Fluent's pseudo-time step (TMR flat
// plate: BiCGStab stalled at its 20-iteration cap, and extra sweeps or
// omega 0.9 diverged). cs bounds the radius by 1
// (Gershgorin over the face coupling); the plate's smoother is then stable
// up to omega 0.9. Only the smoothers use cs; the operator (J x) stays exact.
//
// Packed storage (Metal allows 31 buffer arguments per kernel):
//   fine cell data, 20n floats: [du|dv|dw|cd | Gx|Gy|Gz | Hx|Hy|Hz | Q (9, row major) | cs],
//   G, H and Q unscaled by wc, cs the Schur addition (scaled);
//   coarse levels: a, c per entry, G and H as 3*nnz [x | y | z] (scaled),
//   momentum diagonals and Schur addition as 4n [du | dv | dw | cs].

struct kblock4 { float u, v, w, p; };

KUTTA_FN kblock4 kb4(float u, float v, float w, float p) {
  kblock4 b; b.u = u; b.v = v; b.w = w; b.p = p; return b;
}

// Arrow-shaped diagonal block [[D, e], [h^T, c]], D diagonal. The Schur
// denominator c - h.D^-1.e equals c in interior cells; near boundaries it is
// kept at least half of c so that one cell's block never amplifies.
KUTTA_FN kblock4 coupled_block_solve(kblock4 r, float du, float dv, float dw,
    float ex, float ey, float ez, float hx, float hy, float hz, float c) {
  float iu = 1.0f / du, iv = 1.0f / dv, iw = 1.0f / dw;
  float schur = kmax(c - (hx * ex * iu + hy * ey * iv + hz * ez * iw), 0.5f * c);
  float p = (r.p - hx * iu * r.u - hy * iv * r.v - hz * iw * r.w) / schur;
  return kb4(iu * (r.u - ex * p), iv * (r.v - ey * p), iw * (r.w - ez * p), p);
}

// Neumaier-compensated accumulation; Galerkin diagonals cancel strongly.
struct ksum { float s, e; };
KUTTA_FN ksum ksum_add(ksum a, float b) {
  float t = a.s + b;
  float e = kabs(a.s) >= kabs(b) ? (a.s - t) + b : (b - t) + a.s;
  ksum r; r.s = t; r.e = a.e + e; return r;
}
KUTTA_FN float ksum_value(ksum a) { return a.s + a.e; }
KUTTA_FN ksum ksum0() { ksum r; r.s = 0.0f; r.e = 0.0f; return r; }

// Share of cell c in the face velocity it shares with cell i: the residual's
// momentum-weighted interpolation u_f = (a_i u_i + a_c u_c)/(a_i + a_c)
// (fv_momentum_faces and the face reconstructions).
KUTTA_FN float coupled_share(KIn(float) weight, int i, int c) {
  return weight[c] / kmax(weight[i] + weight[c], 1.0e-30f);
}

// Own-cell continuity-velocity block H_ii (unscaled): the cell's share of
// each internal face velocity plus whole pressure-specified boundary faces,
// whose face velocity is the cell's. Refreshed every solve with the weights.
KUTTA_KERNEL void coupled_velocity_block(int i, KIn(int) offsets, KIn(int) faces, KIn(float) sign,
    KIn(int) adjacent, KIn(int) kind, KIn(float) sx, KIn(float) sy, KIn(float) sz, KIn(float) weight,
    KOut(float) cell, int n, int internal) {
  float hx = 0.0f, hy = 0.0f, hz = 0.0f;
  for (int j = offsets[i]; j < offsets[i + 1]; ++j) {
    int f = faces[j], c = adjacent[j];
    float share = 0.0f;
    if (c >= 0) share = 1.0f - coupled_share(weight, i, c);
    else if (kind[f - internal] == 2 || kind[f - internal] == 4) share = 1.0f;
    float s = share * sign[j];
    hx += s * sx[f]; hy += s * sy[f]; hz += s * sz[f];
  }
  cell[7 * n + i] = hx; cell[8 * n + i] = hy; cell[9 * n + i] = hz;
}

// Volume-flux Rhie-Chow coefficient d(face flux)/d(p_owner - p_other), using
// the same physical diagonal as fv_rhie_chow. Pressure-specified faces
// (kinds 2 and 4) couple to their fixed face value; other boundaries do not.
KUTTA_KERNEL void coupled_face_coefficient(int i, KIn(int) owner, KIn(int) neighbour,
    KIn(int) kind, KIn(float) diffusion, KIn(float) volume, KIn(float) diagonal,
    KOut(float) coefficient, int internal, float scale) {
  int c = owner[i];
  float rc = volume[c] / kmax(diagonal[c], 1.0e-20f);
  if (i < internal) {
    int n = neighbour[i];
    rc = 0.5f * (rc + volume[n] / kmax(diagonal[n], 1.0e-20f));
  } else {
    int b = i - internal;
    if (kind[b] != 2 && kind[b] != 4) { coefficient[i] = 0.0f; return; }
  }
  coefficient[i] = scale * rc * diffusion[i];
}

// Pressure diagonal and RHS. A closed domain doubles the reference cell's
// pressure diagonal: an anchor that lifts the constant pressure null space.
KUTTA_KERNEL void coupled_rhs(int i, KIn(int) offsets, KIn(int) faces, KIn(float) sign,
    KIn(float) coefficient, KIn(float) mass, KIn(float) ru, KIn(float) rv, KIn(float) rw,
    KOut(float) cell, KOut(float) b, int n, int reference, float continuity_scale) {
  ksum c = ksum0();
  float imbalance = 0.0f;
  for (int j = offsets[i]; j < offsets[i + 1]; ++j) {
    int f = faces[j];
    c = ksum_add(c, coefficient[f]);
    imbalance += sign[j] * mass[f];
  }
  float value = ksum_value(c);
  cell[3 * n + i] = i == reference ? 2.0f * value : value;
  b[i] = ru[i]; b[n + i] = rv[i]; b[2 * n + i] = rw[i];
  b[3 * n + i] = -continuity_scale * imbalance;
}

// Jacobian momentum diagonal: the assembled diagonal (pseudo time, no
// implicit relaxation) with the physical part a_P floored at the off-diagonal
// sum, then Fluent's Coupled implicit relaxation a_P / alpha with
// (1 - alpha) / alpha = 1 / CFL (Theory Guide, Steady-State Iterative
// Algorithm, the flow Courant number; eqs. 18.4-23 and 18.4-24 in Fluent 12).
// That is a local time step; the global pseudo time step is negligible in
// cells much smaller than the domain length scale.
//  - The floor: the assembled a_P = sum(D + max(F, 0)) carries the net
//    outflow, so a cell that the current flux only feeds (a sink, common in
//    a non-solenoidal start-up) has a_P near zero, and its row returned
//    updates of tens of U on external-aero meshes. Patankar's rule
//    a_P >= sum |a_nb| keeps the implicit operator an M-matrix.
//  - The Courant term: closed streamlines (separation, a blunt trailing edge's
//    wake, the gap under a wing in ground effect) leave upwind rows with zero
//    row sums, so the momentum block is singular there without it.
// The residual is unchanged and the floor vanishes at convergence (it is the
// continuity residual), so both change only the convergence path.
KUTTA_KERNEL void coupled_copy_diagonal(int i, KIn(int) rows, KIn(int) diag, KIn(float) values,
                                        KIn(float) physical, KOut(float) cell, int offset, float inverse_courant) {
  float off = 0.0f;
  for (int j = rows[i]; j < rows[i + 1]; ++j) if (j != diag[i]) off += kabs(values[j]);
  float floored = kmax(physical[i], off);
  cell[offset + i] = values[diag[i]] + (floored - physical[i]) + inverse_courant * floored;
}

// Absolute row sums of the pressure-force block per component, R_k,c,
// written to r as 3n [x | y | z]. Smoother setup only (see the header).
KUTTA_KERNEL void coupled_fine_gsum(int i, KIn(int) offsets, KIn(int) adjacent, KIn(int) faces,
    KIn(float) sign, KIn(float) sx, KIn(float) sy, KIn(float) sz, KIn(float) dx, KIn(float) dy,
    KIn(float) dz, KIn(float) cell, KOut(float) r, int n, float wc) {
  float q00 = wc * cell[10 * n + i], q01 = wc * cell[11 * n + i], q02 = wc * cell[12 * n + i];
  float q10 = wc * cell[13 * n + i], q11 = wc * cell[14 * n + i], q12 = wc * cell[15 * n + i];
  float q20 = wc * cell[16 * n + i], q21 = wc * cell[17 * n + i], q22 = wc * cell[18 * n + i];
  float half_wc = 0.5f * wc;
  float gx = kabs(wc * cell[4 * n + i]), gy = kabs(wc * cell[5 * n + i]), gz = kabs(wc * cell[6 * n + i]);
  for (int j = offsets[i]; j < offsets[i + 1]; ++j) {
    if (adjacent[j] < 0) continue;
    int f = faces[j];
    float s = sign[j];
    gx += kabs(half_wc * s * sx[f] + q00 * dx[j] + q01 * dy[j] + q02 * dz[j]);
    gy += kabs(half_wc * s * sy[f] + q10 * dx[j] + q11 * dy[j] + q12 * dz[j]);
    gz += kabs(half_wc * s * sz[f] + q20 * dx[j] + q21 * dy[j] + q22 * dz[j]);
  }
  r[i] = gx; r[n + i] = gy; r[2 * n + i] = gz;
}

// cs_i from H_ik = wc s S_f (a_k / (a_i + a_k)) and the neighbours' row sums.
KUTTA_KERNEL void coupled_fine_schur(int i, KIn(int) offsets, KIn(int) adjacent, KIn(int) faces,
    KIn(float) sx, KIn(float) sy, KIn(float) sz, KIn(float) weight, KIn(float) r, KOut(float) cell, int n, float wc) {
  float sum = 0.0f;
  for (int j = offsets[i]; j < offsets[i + 1]; ++j) {
    int c = adjacent[j];
    if (c < 0) continue;
    int f = faces[j];
    sum += wc * coupled_share(weight, i, c) * (kabs(sx[f]) * r[c] / cell[c] + kabs(sy[f]) * r[n + c] / cell[n + c] +
                      kabs(sz[f]) * r[2 * n + c] / cell[2 * n + c]);
  }
  cell[19 * n + i] = sum;
}

// Matrix-free fine-level row: (J x)_i from the face data the momentum and
// Rhie-Chow assemblies already hold, plus the least-squares deltas.
KUTTA_FN kblock4 coupled_fine_row(int i, KIn(int) offsets, KIn(int) faces, KIn(float) sign,
    KIn(int) adjacent, KIn(float) sx, KIn(float) sy, KIn(float) sz, KIn(float) diffusion,
    KIn(float) mu, KIn(float) mass, KIn(float) coefficient, KIn(float) weight, KIn(float) dx, KIn(float) dy,
    KIn(float) dz, KIn(float) cell, KIn(float) x, int n, float wc) {
  float xu = x[i], xv = x[n + i], xw = x[2 * n + i], xp = x[3 * n + i];
  float ex = wc * cell[4 * n + i], ey = wc * cell[5 * n + i], ez = wc * cell[6 * n + i];
  float hx = wc * cell[7 * n + i], hy = wc * cell[8 * n + i], hz = wc * cell[9 * n + i];
  kblock4 y = kb4(cell[i] * xu + ex * xp, cell[n + i] * xv + ey * xp, cell[2 * n + i] * xw + ez * xp,
                  cell[3 * n + i] * xp + hx * xu + hy * xv + hz * xw);
  float q00 = wc * cell[10 * n + i], q01 = wc * cell[11 * n + i], q02 = wc * cell[12 * n + i];
  float q10 = wc * cell[13 * n + i], q11 = wc * cell[14 * n + i], q12 = wc * cell[15 * n + i];
  float q20 = wc * cell[16 * n + i], q21 = wc * cell[17 * n + i], q22 = wc * cell[18 * n + i];
  float half_wc = 0.5f * wc;
  for (int j = offsets[i]; j < offsets[i + 1]; ++j) {
    int c = adjacent[j];
    if (c < 0) continue;
    int f = faces[j];
    float s = sign[j];
    float a = -mu[f] * diffusion[f] + kmin(s * mass[f], 0.0f);
    float bx = half_wc * s * sx[f], by = half_wc * s * sy[f], bz = half_wc * s * sz[f];
    float gx = bx + q00 * dx[j] + q01 * dy[j] + q02 * dz[j];
    float gy = by + q10 * dx[j] + q11 * dy[j] + q12 * dz[j];
    float gz = bz + q20 * dx[j] + q21 * dy[j] + q22 * dz[j];
    float h = wc * s * coupled_share(weight, i, c);
    float pu = x[c], pv = x[n + c], pw = x[2 * n + c], pp = x[3 * n + c];
    y.u += a * pu + gx * pp; y.v += a * pv + gy * pp; y.w += a * pw + gz * pp;
    y.p += h * (sx[f] * pu + sy[f] * pv + sz[f] * pw) - coefficient[f] * pp;
  }
  return y;
}

KUTTA_KERNEL void coupled_fine_spmv(int i, KIn(int) offsets, KIn(int) faces, KIn(float) sign,
    KIn(int) adjacent, KIn(float) sx, KIn(float) sy, KIn(float) sz, KIn(float) diffusion,
    KIn(float) mu, KIn(float) mass, KIn(float) coefficient, KIn(float) weight, KIn(float) dx, KIn(float) dy,
    KIn(float) dz, KIn(float) cell, KIn(float) x, KOut(float) y, int n, float wc) {
  kblock4 r = coupled_fine_row(i, offsets, faces, sign, adjacent, sx, sy, sz, diffusion, mu, mass,
                               coefficient, weight, dx, dy, dz, cell, x, n, wc);
  y[i] = r.u; y[n + i] = r.v; y[2 * n + i] = r.w; y[3 * n + i] = r.p;
}

KUTTA_KERNEL void coupled_fine_residual(int i, KIn(int) offsets, KIn(int) faces, KIn(float) sign,
    KIn(int) adjacent, KIn(float) sx, KIn(float) sy, KIn(float) sz, KIn(float) diffusion,
    KIn(float) mu, KIn(float) mass, KIn(float) coefficient, KIn(float) weight, KIn(float) dx, KIn(float) dy,
    KIn(float) dz, KIn(float) cell, KIn(float) x, KIn(float) b, KOut(float) r, int n, float wc) {
  kblock4 y = coupled_fine_row(i, offsets, faces, sign, adjacent, sx, sy, sz, diffusion, mu, mass,
                               coefficient, weight, dx, dy, dz, cell, x, n, wc);
  r[i] = b[i] - y.u; r[n + i] = b[n + i] - y.v; r[2 * n + i] = b[2 * n + i] - y.w;
  r[3 * n + i] = b[3 * n + i] - y.p;
}

// Damped block Jacobi, double buffered: next = x + omega * B^-1 (b - J x).
// first != 0 means x == 0, so the residual is b and x is not read.
KUTTA_KERNEL void coupled_fine_smooth(int i, KIn(int) offsets, KIn(int) faces, KIn(float) sign,
    KIn(int) adjacent, KIn(float) sx, KIn(float) sy, KIn(float) sz, KIn(float) diffusion,
    KIn(float) mu, KIn(float) mass, KIn(float) coefficient, KIn(float) weight, KIn(float) dx, KIn(float) dy,
    KIn(float) dz, KIn(float) cell, KIn(float) x, KIn(float) b, KOut(float) next,
    int n, float wc, float omega, int first) {
  kblock4 r = kb4(b[i], b[n + i], b[2 * n + i], b[3 * n + i]);
  if (first == 0) {
    kblock4 y = coupled_fine_row(i, offsets, faces, sign, adjacent, sx, sy, sz, diffusion, mu, mass,
                                 coefficient, weight, dx, dy, dz, cell, x, n, wc);
    r = kb4(r.u - y.u, r.v - y.v, r.w - y.w, r.p - y.p);
  }
  kblock4 d = coupled_block_solve(r, cell[i], cell[n + i], cell[2 * n + i],
      wc * cell[4 * n + i], wc * cell[5 * n + i], wc * cell[6 * n + i],
      wc * cell[7 * n + i], wc * cell[8 * n + i], wc * cell[9 * n + i], cell[3 * n + i] + cell[19 * n + i]);
  float base_u = first == 0 ? x[i] : 0.0f, base_v = first == 0 ? x[n + i] : 0.0f;
  float base_w = first == 0 ? x[2 * n + i] : 0.0f, base_p = first == 0 ? x[3 * n + i] : 0.0f;
  next[i] = base_u + omega * d.u; next[n + i] = base_v + omega * d.v;
  next[2 * n + i] = base_w + omega * d.w; next[3 * n + i] = base_p + omega * d.p;
}

// Coarse levels: CSR with a (uu), c (pp), G (up) and H (pu) per entry and the
// three momentum diagonals per row. G and H already include wc.
KUTTA_FN kblock4 coupled_coarse_row(int i, KIn(int) rows, KIn(int) cols, KIn(int) diag,
    KIn(float) a, KIn(float) c, KIn(float) g, KIn(float) h, KIn(float) d, KIn(float) x,
    int n, int nnz) {
  int k = diag[i];
  float xu = x[i], xv = x[n + i], xw = x[2 * n + i], xp = x[3 * n + i];
  kblock4 y = kb4(d[i] * xu + g[k] * xp, d[n + i] * xv + g[nnz + k] * xp,
                  d[2 * n + i] * xw + g[2 * nnz + k] * xp,
                  c[k] * xp + h[k] * xu + h[nnz + k] * xv + h[2 * nnz + k] * xw);
  for (int j = rows[i]; j < rows[i + 1]; ++j) {
    if (j == k) continue;
    int m = cols[j];
    float pu = x[m], pv = x[n + m], pw = x[2 * n + m], pp = x[3 * n + m];
    y.u += a[j] * pu + g[j] * pp; y.v += a[j] * pv + g[nnz + j] * pp; y.w += a[j] * pw + g[2 * nnz + j] * pp;
    y.p += h[j] * pu + h[nnz + j] * pv + h[2 * nnz + j] * pw + c[j] * pp;
  }
  return y;
}

KUTTA_KERNEL void coupled_coarse_residual(int i, KIn(int) rows, KIn(int) cols, KIn(int) diag,
    KIn(float) a, KIn(float) c, KIn(float) g, KIn(float) h, KIn(float) d, KIn(float) x,
    KIn(float) b, KOut(float) r, int n, int nnz) {
  kblock4 y = coupled_coarse_row(i, rows, cols, diag, a, c, g, h, d, x, n, nnz);
  r[i] = b[i] - y.u; r[n + i] = b[n + i] - y.v; r[2 * n + i] = b[2 * n + i] - y.w;
  r[3 * n + i] = b[3 * n + i] - y.p;
}

KUTTA_KERNEL void coupled_coarse_smooth(int i, KIn(int) rows, KIn(int) cols, KIn(int) diag,
    KIn(float) a, KIn(float) c, KIn(float) g, KIn(float) h, KIn(float) d, KIn(float) x,
    KIn(float) b, KOut(float) next, int n, int nnz, float omega, int first) {
  kblock4 r = kb4(b[i], b[n + i], b[2 * n + i], b[3 * n + i]);
  if (first == 0) {
    kblock4 y = coupled_coarse_row(i, rows, cols, diag, a, c, g, h, d, x, n, nnz);
    r = kb4(r.u - y.u, r.v - y.v, r.w - y.w, r.p - y.p);
  }
  int k = diag[i];
  kblock4 s = coupled_block_solve(r, d[i], d[n + i], d[2 * n + i], g[k], g[nnz + k], g[2 * nnz + k],
                                  h[k], h[nnz + k], h[2 * nnz + k], c[k] + d[3 * n + i]);
  float base_u = first == 0 ? x[i] : 0.0f, base_v = first == 0 ? x[n + i] : 0.0f;
  float base_w = first == 0 ? x[2 * n + i] : 0.0f, base_p = first == 0 ? x[3 * n + i] : 0.0f;
  next[i] = base_u + omega * s.u; next[n + i] = base_v + omega * s.v;
  next[2 * n + i] = base_w + omega * s.w; next[3 * n + i] = base_p + omega * s.p;
}

// Linear search in one sorted coarse row; rows hold a few dozen entries.
KUTTA_FN int coupled_entry(KIn(int) rows, KIn(int) cols, int row, int col) {
  for (int j = rows[row]; j < rows[row + 1]; ++j) if (cols[j] == col) return j;
  return -1;
}

// Piecewise-constant Galerkin from the matrix-free fine level. One thread owns
// one coarse row and all of its entries, so the sum order is fixed.
KUTTA_KERNEL void coupled_galerkin_fine(int i, KIn(int) member_offsets, KIn(int) members,
    KIn(int) aggregate, KIn(int) offsets, KIn(int) faces, KIn(float) sign, KIn(int) adjacent,
    KIn(float) sx, KIn(float) sy, KIn(float) sz, KIn(float) diffusion, KIn(float) mu,
    KIn(float) mass, KIn(float) coefficient, KIn(float) weight, KIn(float) dx, KIn(float) dy, KIn(float) dz,
    KIn(float) cell, KIn(int) rows, KIn(int) cols,
    KOut(float) a, KOut(float) c, KOut(float) g, KOut(float) h, KOut(float) d,
    int nf, int n, int nnz, float wc) {
  for (int j = rows[i]; j < rows[i + 1]; ++j) {
    a[j] = 0.0f; c[j] = 0.0f;
    g[j] = 0.0f; g[nnz + j] = 0.0f; g[2 * nnz + j] = 0.0f;
    h[j] = 0.0f; h[nnz + j] = 0.0f; h[2 * nnz + j] = 0.0f;
  }
  ksum su = ksum0(), sv = ksum0(), sw = ksum0(), sc = ksum0();
  ksum sgx = ksum0(), sgy = ksum0(), sgz = ksum0(), shx = ksum0(), shy = ksum0(), shz = ksum0();
  float half_wc = 0.5f * wc;
  for (int m = member_offsets[i]; m < member_offsets[i + 1]; ++m) {
    int q = members[m];
    su = ksum_add(su, cell[q]); sv = ksum_add(sv, cell[nf + q]); sw = ksum_add(sw, cell[2 * nf + q]);
    sc = ksum_add(sc, cell[3 * nf + q]);
    sgx = ksum_add(sgx, wc * cell[4 * nf + q]); sgy = ksum_add(sgy, wc * cell[5 * nf + q]);
    sgz = ksum_add(sgz, wc * cell[6 * nf + q]);
    shx = ksum_add(shx, wc * cell[7 * nf + q]); shy = ksum_add(shy, wc * cell[8 * nf + q]);
    shz = ksum_add(shz, wc * cell[9 * nf + q]);
    float q00 = wc * cell[10 * nf + q], q01 = wc * cell[11 * nf + q], q02 = wc * cell[12 * nf + q];
    float q10 = wc * cell[13 * nf + q], q11 = wc * cell[14 * nf + q], q12 = wc * cell[15 * nf + q];
    float q20 = wc * cell[16 * nf + q], q21 = wc * cell[17 * nf + q], q22 = wc * cell[18 * nf + q];
    for (int j = offsets[q]; j < offsets[q + 1]; ++j) {
      int other = adjacent[j];
      if (other < 0) continue;
      int f = faces[j];
      float s = sign[j];
      float uu = -mu[f] * diffusion[f] + kmin(s * mass[f], 0.0f);
      float bx = half_wc * s * sx[f], by = half_wc * s * sy[f], bz = half_wc * s * sz[f];
      float gx = bx + q00 * dx[j] + q01 * dy[j] + q02 * dz[j];
      float gy = by + q10 * dx[j] + q11 * dy[j] + q12 * dz[j];
      float gz = bz + q20 * dx[j] + q21 * dy[j] + q22 * dz[j];
      float pp = -coefficient[f];
      float share = wc * s * coupled_share(weight, q, other);
      float hx = share * sx[f], hy = share * sy[f], hz = share * sz[f];
      int target = aggregate[other];
      if (target == i) {
        su = ksum_add(su, uu); sv = ksum_add(sv, uu); sw = ksum_add(sw, uu);
        sgx = ksum_add(sgx, gx); sgy = ksum_add(sgy, gy); sgz = ksum_add(sgz, gz);
        shx = ksum_add(shx, hx); shy = ksum_add(shy, hy); shz = ksum_add(shz, hz);
        sc = ksum_add(sc, pp);
      } else {
        int e = coupled_entry(rows, cols, i, target);
        a[e] += uu; c[e] += pp;
        g[e] += gx; g[nnz + e] += gy; g[2 * nnz + e] += gz;
        h[e] += hx; h[nnz + e] += hy; h[2 * nnz + e] += hz;
      }
    }
  }
  // Metal's 31-argument limit leaves no slot for the diagonal index array.
  int k = coupled_entry(rows, cols, i, i);
  d[i] = ksum_value(su); d[n + i] = ksum_value(sv); d[2 * n + i] = ksum_value(sw);
  c[k] = ksum_value(sc);
  g[k] = ksum_value(sgx); g[nnz + k] = ksum_value(sgy); g[2 * nnz + k] = ksum_value(sgz);
  h[k] = ksum_value(shx); h[nnz + k] = ksum_value(shy); h[2 * nnz + k] = ksum_value(shz);
}

KUTTA_KERNEL void coupled_galerkin_coarse(int i, KIn(int) member_offsets, KIn(int) members,
    KIn(int) aggregate, KIn(int) frows, KIn(int) fcols, KIn(int) fdiag,
    KIn(float) fa, KIn(float) fc, KIn(float) fg, KIn(float) fh, KIn(float) fd,
    KIn(int) rows, KIn(int) cols, KIn(int) diag,
    KOut(float) a, KOut(float) c, KOut(float) g, KOut(float) h, KOut(float) d,
    int nf, int fnnz, int n, int nnz) {
  for (int j = rows[i]; j < rows[i + 1]; ++j) {
    a[j] = 0.0f; c[j] = 0.0f;
    g[j] = 0.0f; g[nnz + j] = 0.0f; g[2 * nnz + j] = 0.0f;
    h[j] = 0.0f; h[nnz + j] = 0.0f; h[2 * nnz + j] = 0.0f;
  }
  ksum su = ksum0(), sv = ksum0(), sw = ksum0(), sc = ksum0();
  ksum sgx = ksum0(), sgy = ksum0(), sgz = ksum0(), shx = ksum0(), shy = ksum0(), shz = ksum0();
  for (int m = member_offsets[i]; m < member_offsets[i + 1]; ++m) {
    int row = members[m];
    su = ksum_add(su, fd[row]); sv = ksum_add(sv, fd[nf + row]); sw = ksum_add(sw, fd[2 * nf + row]);
    for (int j = frows[row]; j < frows[row + 1]; ++j) {
      int target = aggregate[fcols[j]];
      if (j == fdiag[row] || target == i) {
        if (j != fdiag[row]) { su = ksum_add(su, fa[j]); sv = ksum_add(sv, fa[j]); sw = ksum_add(sw, fa[j]); }
        sgx = ksum_add(sgx, fg[j]); sgy = ksum_add(sgy, fg[fnnz + j]); sgz = ksum_add(sgz, fg[2 * fnnz + j]);
        shx = ksum_add(shx, fh[j]); shy = ksum_add(shy, fh[fnnz + j]); shz = ksum_add(shz, fh[2 * fnnz + j]);
        sc = ksum_add(sc, fc[j]);
      } else {
        int e = coupled_entry(rows, cols, i, target);
        a[e] += fa[j]; c[e] += fc[j];
        g[e] += fg[j]; g[nnz + e] += fg[fnnz + j]; g[2 * nnz + e] += fg[2 * fnnz + j];
        h[e] += fh[j]; h[nnz + e] += fh[fnnz + j]; h[2 * nnz + e] += fh[2 * fnnz + j];
      }
    }
  }
  int k = diag[i];
  d[i] = ksum_value(su); d[n + i] = ksum_value(sv); d[2 * n + i] = ksum_value(sw);
  c[k] = ksum_value(sc);
  g[k] = ksum_value(sgx); g[nnz + k] = ksum_value(sgy); g[2 * nnz + k] = ksum_value(sgz);
  h[k] = ksum_value(shx); h[nnz + k] = ksum_value(shy); h[2 * nnz + k] = ksum_value(shz);
}

// Coarse-level row sums R_j,c of |G| (diagonal included), as 3n in r.
KUTTA_KERNEL void coupled_coarse_gsum(int i, KIn(int) rows, KIn(float) g, KOut(float) r, int n, int nnz) {
  float gx = 0.0f, gy = 0.0f, gz = 0.0f;
  for (int j = rows[i]; j < rows[i + 1]; ++j) {
    gx += kabs(g[j]); gy += kabs(g[nnz + j]); gz += kabs(g[2 * nnz + j]);
  }
  r[i] = gx; r[n + i] = gy; r[2 * n + i] = gz;
}

// Coarse-level cs_i from the Galerkin H entries. Aggregated momentum
// diagonals stay positive (pseudo time); the floor only guards round-off.
KUTTA_KERNEL void coupled_coarse_schur(int i, KIn(int) rows, KIn(int) cols, KIn(int) diag,
    KIn(float) h, KIn(float) r, KOut(float) d, int n, int nnz) {
  float sum = 0.0f;
  for (int j = rows[i]; j < rows[i + 1]; ++j) {
    if (j == diag[i]) continue;
    int m = cols[j];
    sum += kabs(h[j]) * r[m] / kmax(d[m], 1.0e-30f) + kabs(h[nnz + j]) * r[n + m] / kmax(d[n + m], 1.0e-30f) +
           kabs(h[2 * nnz + j]) * r[2 * n + m] / kmax(d[2 * n + m], 1.0e-30f);
  }
  d[3 * n + i] = sum;
}

KUTTA_KERNEL void coupled_restrict(int i, KIn(int) member_offsets, KIn(int) members,
    KIn(float) fine, KOut(float) coarse, int nf, int nc) {
  for (int k = 0; k < 4; ++k) {
    float sum = 0.0f;
    for (int m = member_offsets[i]; m < member_offsets[i + 1]; ++m) sum += fine[k * nf + members[m]];
    coarse[k * nc + i] = sum;
  }
}

KUTTA_KERNEL void coupled_prolong(int i, KIn(int) aggregate, KIn(float) coarse, KOut(float) fine,
    int nf, int nc) {
  int g = aggregate[i];
  for (int k = 0; k < 4; ++k) fine[k * nf + i] += coarse[k * nc + g];
}

// Fluent's explicit relaxation of the coupled update (Solution Controls,
// pseudo time explicit relaxation factors), then the closed-domain gauge.
KUTTA_KERNEL void coupled_update(int i, KIn(float) delta, KOut(float) u, KOut(float) v,
    KOut(float) w, KOut(float) p, int n, float velocity_relaxation, float pressure_relaxation) {
  u[i] += velocity_relaxation * delta[i];
  v[i] += velocity_relaxation * delta[n + i];
  w[i] += velocity_relaxation * delta[2 * n + i];
  p[i] += pressure_relaxation * delta[3 * n + i];
}

KUTTA_KERNEL void coupled_gauge(int i, KIn(float) reference_value, KOut(float) p) {
  p[i] -= reference_value[0];
}
KUTTA_KERNEL void coupled_pick(int i, KIn(float) p, KOut(float) out, int index) {
  if (i == 0) out[0] = p[index];
}
