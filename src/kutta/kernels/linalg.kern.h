#include "kutta/kernels/prelude.h"

// One thread owns a row. Long unstructured rows would benefit from a
// warp/SIMD-per-row variant; this gather baseline also handles short FV rows.
KUTTA_KERNEL void csr_spmv(int i, KIn(int) rows, KIn(int) cols, KIn(float) values,
                          KIn(float) x, KOut(float) y) {
  float sum = 0.0f;
  float error = 0.0f;
  for (int j = rows[i]; j < rows[i + 1]; ++j) {
    float term = values[j] * x[cols[j]] - error;
    float next = sum + term;
    error = (next - sum) - term;
    sum = next;
  }
  y[i] = sum;
}

// Dekker product correction uses only float operations (also on Metal). This
// matters when a tiny residual subtracts nearly equal stencil products. Fall
// back to a rounded product outside the safe splitting range.
KUTTA_FN float csr_product_error(float a, float b, float product) {
  if (kabs(a) > 1.0e30f || kabs(b) > 1.0e30f) return 0.0f;
  float ca = 4097.0f * a;
  float cb = 4097.0f * b;
  float ah = ca - (ca - a);
  float bh = cb - (cb - b);
  float al = a - ah;
  float bl = b - bh;
  return ((ah * bh - product) + ah * bl + al * bh) + al * bl;
}

KUTTA_FN float csr_row_residual(int i, KIn(int) rows, KIn(int) cols, KIn(float) values,
                                 KIn(float) x, KIn(float) b) {
  float sum = b[i];
  float error = 0.0f;
  for (int j = rows[i]; j < rows[i + 1]; ++j) {
    float product = values[j] * x[cols[j]];
    float term = -product;
    float next = sum + term;
    error += kabs(sum) >= kabs(term) ? (sum - next) + term : (term - next) + sum;
    error -= csr_product_error(values[j], x[cols[j]], product);
    sum = next;
  }
  return sum + error;
}

KUTTA_KERNEL void csr_residual(int i, KIn(int) rows, KIn(int) cols, KIn(float) values,
                              KIn(float) x, KIn(float) b, KOut(float) r) {
  r[i] = csr_row_residual(i, rows, cols, values, x, b);
}

KUTTA_KERNEL void krylov_diagonal(int i, KIn(int) rows, KIn(int) diag,
                                 KIn(float) values, KOut(float) inverse,
                                 KOut(float) invalid, int use_l1) {
  float d = values[diag[i]];
  float off = 0.0f;
  float bad = 0.0f;
  for (int j = rows[i]; j < rows[i + 1]; ++j) {
    float v = values[j];
    if (v - v != 0.0f) bad = 1.0f;
    if (j != diag[i]) off += kabs(v);
  }
  float inv = d != 0.0f ? 1.0f / d : 0.0f;
  // D_l1 >= (|D| + sum |offdiag|)/2 bounds the spectrum by two for SPD
  // operators. Symmetric pre/post sweeps therefore give an SPD V-cycle.
  float l1 = kmax(kabs(d), 0.5f * (kabs(d) + off));
  inverse[i] = use_l1 != 0 ? (l1 != 0.0f ? 1.0f / l1 : 0.0f) : inv;
  invalid[i] = (bad != 0.0f || inv - inv != 0.0f || l1 - l1 != 0.0f) ?
      2.0f : (use_l1 != 0 && d < 0.0f ? 3.0f : (d == 0.0f ? 1.0f : 0.0f));
}

KUTTA_KERNEL void krylov_multiply(int i, KIn(float) r, KIn(float) inverse,
                                 KOut(float) z) {
  z[i] = r[i] * inverse[i];
}

KUTTA_KERNEL void krylov_smooth(int i, KIn(float) r, KIn(float) inverse,
                               KOut(float) x, float omega) {
  x[i] += omega * inverse[i] * r[i];
}

KUTTA_KERNEL void krylov_shift(int i, KOut(float) x, float mean) { x[i] -= mean; }

KUTTA_KERNEL void krylov_direction(int i, KIn(float) z, KOut(float) p, float beta) {
  p[i] = z[i] + beta * p[i];
}

KUTTA_KERNEL void krylov_step(int i, KOut(float) x, KOut(float) r,
                             KIn(float) p, KIn(float) ap, float alpha) {
  x[i] += alpha * p[i];
  r[i] -= alpha * ap[i];
}

KUTTA_KERNEL void krylov_bicg_direction(int i, KIn(float) r, KIn(float) v,
                                       KOut(float) p, float beta, float omega) {
  p[i] = r[i] + beta * (p[i] - omega * v[i]);
}

KUTTA_KERNEL void krylov_difference(int i, KIn(float) a, KIn(float) b,
                                   KOut(float) out, float scale) {
  out[i] = a[i] - scale * b[i];
}

KUTTA_KERNEL void krylov_bicg_step(int i, KOut(float) x, KIn(float) p,
                                  KIn(float) s, float alpha, float omega) {
  x[i] += alpha * p[i] + omega * s[i];
}

// Frozen Krylov batches still enqueue dispatches; guards prevent their large
// memory sweeps after a convergence candidate or breakdown has been detected.
KUTTA_KERNEL void csr_spmv_guarded(int i, KIn(int) rows, KIn(int) cols, KIn(float) values,
                                  KIn(float) x, KOut(float) y, KIn(int) status) {
  if (status[0] == 0) csr_spmv(i, rows, cols, values, x, y);
}
KUTTA_KERNEL void csr_residual_guarded(int i, KIn(int) rows, KIn(int) cols, KIn(float) values,
                                      KIn(float) x, KIn(float) b, KOut(float) r, KIn(int) status) {
  if (status[0] == 0) csr_residual(i, rows, cols, values, x, b, r);
}
KUTTA_KERNEL void krylov_multiply_guarded(int i, KIn(float) r, KIn(float) inverse,
                                         KOut(float) z, KIn(int) status) {
  if (status[0] == 0) krylov_multiply(i, r, inverse, z);
}
// One owner for an entire <=24-row coarse system. Local arrays preserve the
// simultaneous Jacobi updates and exactly the ordinary residual's arithmetic.
// This replaces 64 tiny dispatches, with no inter-thread barrier or host solve.
KUTTA_KERNEL void amg_coarse_solve(int i, KIn(int) rows, KIn(int) cols,
                                   KIn(float) values, KOut(float) x, KIn(float) b,
                                   KIn(float) inverse, KIn(int) status, int count, int sweeps) {
  if (i != 0 || status[0] != 0) return;
  float current[24], r[24];
  for (int row = 0; row < count; ++row) current[row] = 0.0f;
  for (int sweep = 0; sweep < sweeps; ++sweep) {
    for (int row = 0; row < count; ++row) {
      float sum = b[row], error = 0.0f;
      for (int j = rows[row]; j < rows[row + 1]; ++j) {
        float product = values[j] * current[cols[j]], term = -product;
        float next = sum + term;
        error += kabs(sum) >= kabs(term) ? (sum - next) + term : (term - next) + sum;
        error -= csr_product_error(values[j], current[cols[j]], product);
        sum = next;
      }
      r[row] = sum + error;
    }
    for (int row = 0; row < count; ++row) current[row] += (2.0f / 3.0f) * inverse[row] * r[row];
  }
  for (int row = 0; row < count; ++row) x[row] = current[row];
}


KUTTA_KERNEL void amg_smooth(int i, KIn(int) rows, KIn(int) cols, KIn(float) values,
                             KIn(float) x, KIn(float) b, KIn(float) inverse,
                             KOut(float) next, KIn(int) status) {
  if (status[0] != 0) return;
  float r = csr_row_residual(i, rows, cols, values, x, b);
  next[i] = x[i] + (2.0f / 3.0f) * inverse[i] * r;
}
