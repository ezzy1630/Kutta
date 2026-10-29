#pragma once
#include "kutta/kernels/prelude.h"

// Fluent stores sigma as a divisor, not the multiplier used by NASA TMR.
KUTTA_FN float sst_inv_sigma_k(float f1) {
  return f1 / 1.176f + (1.0f - f1);
}
KUTTA_FN float sst_inv_sigma_w(float f1) {
  return f1 / 2.0f + (1.0f - f1) / 1.168f;
}
KUTTA_FN float sst_beta(float f1) { return 0.075f * f1 + 0.0828f * (1.0f - f1); }
KUTTA_FN float sst_gamma(float f1, int set_2003) {
  float g1 = set_2003 ? 5.0f / 9.0f : 0.075f / 0.09f - 0.41f * 0.41f / (2.0f * 0.3f);
  float g2 = set_2003 ? 0.44f : 0.0828f / 0.09f - 0.41f * 0.41f / (1.168f * 0.3f);
  return f1 * g1 + (1.0f - f1) * g2;
}
KUTTA_FN float sst_alpha_star(float ret, int low_re) {
  // TG 4-75..79: alpha*_0 = beta_i/3 from STANDARD k-omega (0.072).
  return low_re ? 1.0f - 0.976f / (1.0f + ret / 6.0f) : 1.0f;
}
KUTTA_FN float sst_beta_star(float ret, int low_re) {
  if (!low_re) return 0.09f;
  // Ratio form avoids overflow in Re_T^4 at the high-Re limit (TG 4-93).
  float r = kmin(ret / 8.0f, 100.0f);
  float r4 = r * r * r * r;
  return 0.09f * ((4.0f / 15.0f) + r4) / (1.0f + r4);
}
KUTTA_FN float sst_power_blend(float a, float b, float p) {
  float hi = kmax(a, b), lo = kmin(a, b);
  return hi == 0.0f ? 0.0f : hi * kpow(1.0f + kpow(lo / hi, p), 1.0f / p);
}

KUTTA_KERNEL void sst_rates(int i,
    KIn(float) ux, KIn(float) uy, KIn(float) uz,
    KIn(float) vx, KIn(float) vy, KIn(float) vz,
    KIn(float) wx, KIn(float) wy, KIn(float) wz,
    KOut(float) strain, KOut(float) vorticity) {
  float xy = uy[i] + vx[i], xz = uz[i] + wx[i], yz = vz[i] + wy[i];
  strain[i] = ksqrt(2.0f * (ux[i]*ux[i] + vy[i]*vy[i] + wz[i]*wz[i])
                    + xy*xy + xz*xz + yz*yz);
  xy = uy[i] - vx[i]; xz = uz[i] - wx[i]; yz = vz[i] - wy[i];
  vorticity[i] = ksqrt(xy*xy + xz*xz + yz*yz);
}

KUTTA_KERNEL void sst_blending(int i,
    KIn(float) rho, KIn(float) mu, KIn(float) k, KIn(float) omega, KIn(float) distance,
    KIn(float) kx, KIn(float) ky, KIn(float) kz,
    KIn(float) ox, KIn(float) oy, KIn(float) oz,
    KOut(float) f1, KOut(float) f2, KOut(float) cdkw,
    float cd_floor, float k_floor, float omega_floor) {
  float ke = kmax(k[i], k_floor), w = kmax(omega[i], omega_floor), y = distance[i];
  float cd = 2.0f * rho[i] * (kx[i]*ox[i] + ky[i]*oy[i] + kz[i]*oz[i]) / (1.168f * w);
  float a = ksqrt(ke) / (0.09f * w * y);
  float b = 500.0f * (mu[i] / rho[i]) / (y*y*w);
  float c = 4.0f * rho[i] * ke / (1.168f * kmax(cd, cd_floor) * y*y);
  // Saturating tanh before powers prevents overflow without changing float F1/F2.
  float p1 = kmin(kmin(kmax(a, b), c), 10.0f);
  float p2 = kmin(kmax(2.0f*a, b), 10.0f);
  f1[i] = ktanh(p1*p1*p1*p1);
  f2[i] = ktanh(p2*p2);
  cdkw[i] = cd;  // negative cross diffusion is a real omega sink
}

KUTTA_KERNEL void sst_viscosity(int i,
    KIn(float) rho, KIn(float) mu, KIn(float) k, KIn(float) omega,
    KIn(float) strain, KIn(float) vorticity, KIn(float) f1, KIn(float) f2,
    KOut(float) mut, KOut(float) gamma_k, KOut(float) gamma_w,
    int low_re, int use_vorticity, float max_ratio, float k_floor, float omega_floor) {
  float ke = kmax(k[i], k_floor), w = kmax(omega[i], omega_floor);
  float ret = (rho[i] / mu[i]) * (ke / w);
  float alpha = sst_alpha_star(ret, low_re);
  float rate = use_vorticity ? vorticity[i] : strain[i];
  float mt = rho[i] * ke / kmax(w / alpha, rate * f2[i] / 0.31f);
  mt = kmin(mt, max_ratio * mu[i]);
  mut[i] = mt;
  gamma_k[i] = mu[i] + mt * sst_inv_sigma_k(f1[i]);
  gamma_w[i] = mu[i] + mt * sst_inv_sigma_w(f1[i]);
}

KUTTA_KERNEL void sst_wall(int i,
    KIn(int) owner, KIn(int) wall, int n_internal,
    KIn(float) ax, KIn(float) ay, KIn(float) az,
    KIn(float) ux, KIn(float) uy, KIn(float) uz,
    KIn(float) wall_ux, KIn(float) wall_uy, KIn(float) wall_uz,
    KIn(float) dy, KIn(float) rho, KIn(float) mu, KIn(float) k, KIn(float) omega,
    KOut(float) tau, KOut(float) shear_x, KOut(float) shear_y, KOut(float) shear_z,
    KOut(float) u_star, KOut(float) u_tau, KOut(float) y_plus, KOut(float) y_plus_utau,
    KOut(float) fixed_omega, KOut(float) production,
    float k_floor, float omega_floor) {
  tau[i] = 0.0f; shear_x[i] = 0.0f; shear_y[i] = 0.0f; shear_z[i] = 0.0f;
  u_star[i] = 0.0f; u_tau[i] = 0.0f; y_plus[i] = 0.0f; y_plus_utau[i] = 0.0f;
  fixed_omega[i] = 0.0f; production[i] = 0.0f;
  if (!wall[i]) return;
  int f = n_internal + i, c = owner[f];
  float r = rho[c], m = mu[c], d = dy[i], ke = kmax(k[c], k_floor);
  float w = kmax(omega[c], omega_floor);
  kvec3 normal = kv3(ax[f], ay[f], az[f]);
  normal = (1.0f / knorm(normal)) * normal;
  kvec3 relative = kv3(ux[c]-wall_ux[i], uy[c]-wall_uy[i], uz[c]-wall_uz[i]);
  kvec3 tangent = relative - kdot(relative, normal) * normal;
  float speed = knorm(tangent);
  float us = ksqrt((m / r) * speed / d + 0.3f * ke);
  float yp = (r / m) * us * d;
  float log_up = klog(9.793f * kmax(yp, 0.2f)) / 0.4187f;
  float ratio = kmin(yp, log_up) / kmax(yp, log_up);
  float blend = kpow(1.0f + ratio*ratio*ratio*ratio, 0.25f);
  // Equivalent to rho*u*U/u+, but defined at U=k=0 and at y+ -> 0.
  float stress = (yp <= log_up ? m * speed / d : r * us * speed / log_up) * blend;
  float ut = ksqrt(stress / r);
  tau[i] = stress; u_star[i] = us; u_tau[i] = ut;
  y_plus[i] = yp; y_plus_utau[i] = (r / m) * ut * d;
  if (speed > 0.0f) {
    kvec3 shear = (-stress / speed) * tangent;
    shear_x[i] = shear.x; shear_y[i] = shear.y; shear_z[i] = shear.z;
  }
  // Dimensional branches avoid the removable u*^2/y+^2 singularity.
  float w_lam = 2.0f * (m / r) / (0.072f * d*d);
  float w_log = us / (0.3f * 0.4187f * d);
  fixed_omega[i] = sst_power_blend(w_lam, w_log, 1.3f);
  float p_lam = (r * ke / w) * (stress / m) * (stress / m);
  // Analytic log-law derivative 1/(kappa*y+), TG 4-421; clamp in u+ is
  // for its VALUE (4-416). The y+ -> 0 production limit is zero.
  float p_log = yp > 0.0f ? (stress*stress / m) / (0.4187f * yp) : 0.0f;
  float hi = kmax(p_lam, p_log), lo = kmin(p_lam, p_log);
  production[i] = hi > 0.0f ? lo / (1.0f + lo / hi) : 0.0f;
}

KUTTA_KERNEL void sst_wall_cells(int i,
    KIn(int) offsets, KIn(int) faces, KIn(int) wall, KIn(float) dy,
    KIn(float) face_omega, KIn(float) face_production, int n_internal,
    KOut(int) fixed_mask, KOut(int) selected_face,
    KOut(float) fixed_omega, KOut(float) production) {
  int selected = -1;
  // Multi-wall cells are unresolved in Fluent's reference case. Use the
  // nearest wall plane, with an explicit deterministic tie break.
  for (int j = offsets[i]; j < offsets[i+1]; ++j) {
    int b = faces[j] - n_internal;
    if (b >= 0 && wall[b]) {
      if (selected < 0 || dy[b] < dy[selected] || (dy[b] == dy[selected] && b < selected))
        selected = b;
    }
  }
  selected_face[i] = selected;
  fixed_mask[i] = selected >= 0 ? 1 : 0;
  fixed_omega[i] = selected >= 0 ? face_omega[selected] : 0.0f;
  production[i] = selected >= 0 ? face_production[selected] : 0.0f;
}

KUTTA_KERNEL void sst_sources(int i,
    KIn(float) rho, KIn(float) mu, KIn(float) k, KIn(float) omega,
    KIn(float) strain, KIn(float) vorticity, KIn(float) f1, KIn(float) cdkw,
    KIn(float) mut, KIn(int) wall_mask, KIn(float) wall_pk,
    KOut(float) k_source, KOut(float) k_sink,
    KOut(float) omega_source, KOut(float) omega_sink, KOut(float) production_k,
    int low_re, int gamma_2003, int production_mode, int limiter,
    int limit_omega, int wall_override, float clip, float k_floor, float omega_floor) {
  float r = rho[i], ke = kmax(k[i], k_floor), w = kmax(omega[i], omega_floor);
  float ret = (r / mu[i]) * (ke / w);
  float beta_star = sst_beta_star(ret, low_re);
  float rate2 = production_mode == 1 ? strain[i]*vorticity[i] :
                (production_mode == 2 ? vorticity[i]*vorticity[i] : strain[i]*strain[i]);
  float raw_pk = mut[i] * rate2;
  if (wall_override && wall_mask[i]) raw_pk = wall_pk[i];
  float cap = clip * r * beta_star * ke * w;
  float pk = limiter ? kmin(raw_pk, cap) : raw_pk;
  float pw_k = limiter && limit_omega ? pk : raw_pk;
  float gamma = sst_gamma(f1[i], gamma_2003);
  // TG 4-83/84, 4-110/111: alpha* cancels in alpha*alpha; do not
  // damp G_omega twice when low-Re corrections are enabled.
  if (low_re) gamma *= 1.0f - (8.0f/9.0f) / (1.0f + ret/2.95f);
  float pw = mut[i] > 0.0f ? gamma * r * (pw_k / mut[i]) : gamma * r * rate2;
  float cross = (1.0f - f1[i]) * cdkw[i];
  k_source[i] = pk;
  k_sink[i] = r * beta_star * w;
  // Patankar linearization: negative CDkw becomes an implicit sink;
  // positive CDkw stays on the RHS. Both coefficients remain nonnegative.
  omega_source[i] = pw + kmax(cross, 0.0f);
  omega_sink[i] = r * sst_beta(f1[i]) * w + kmax(-cross, 0.0f) / w;
  production_k[i] = pk;
}

KUTTA_KERNEL void sst_face_diffusion(int i,
    KIn(int) owner, KIn(int) neighbour, KIn(float) weight, KIn(int) wall,
    KIn(float) mu, KIn(float) gamma_k, KIn(float) gamma_w, int n_internal,
    KOut(float) face_k, KOut(float) face_w) {
  int c = owner[i];
  if (i < n_internal) {
    int n = neighbour[i]; float a = weight[i];
    face_k[i] = a * gamma_k[c] + (1.0f - a) * gamma_k[n];
    face_w[i] = a * gamma_w[c] + (1.0f - a) * gamma_w[n];
  } else {
    // mu_t vanishes on a smooth wall. The assembler still enforces zero k
    // flux and fixes wall-cell omega; other boundary coefficients extrapolate.
    face_k[i] = wall[i-n_internal] ? mu[c] : gamma_k[c];
    face_w[i] = wall[i-n_internal] ? mu[c] : gamma_w[c];
  }
}
