#pragma once
#include "kutta/kernels/prelude.h"
#include "kutta/kernels/sst.kern.h"

// Fluent's Intermittency Transition Model (one-equation gamma, Menter et al.
// 2015) coupled to SST. Equation numbers are Fluent 2026 R1 Theory Guide
// section 4.8.2 (gamma) and 4-205..207 (k coupling, shared with Transition
// SST). Terms Fluent names but does not print (P_k^lim, F1 = max(F1, F3)) are
// Menter et al. 2015 as restated by NASA TMR's "SST-2003-Menter-gamma-2015"
// page. src/kutta/fv/turbulence/transition.md lists every source.

// Fluent's rescaled pressure-gradient parameter, Eq. 4-216 and its bound 4-217.
// It is 14.68x the paper's lambda, so C_PG1 = 1 and C_PG2 = -.5 are the paper's
// 14.68 and -7.34 in these units.
KUTTA_FN float transition_lambda(float dvdy, float d, float nu) {
  float lambda = -0.1111f * dvdy * (d * d / nu) + 0.1875f;
  return kmin(kmax(lambda, -10.0f), 10.0f);
}
// F_PG with the paper's robustness floor max(F_PG, 0); with Fluent's default
// C_PG3 = 0 the floor never binds.
KUTTA_FN float transition_fpg(float lambda, float cpg1, float cpg2, float cpg3,
                              float cpg1_lim, float cpg2_lim) {
  float f = lambda >= 0.0f ? kmin(1.0f + cpg1 * lambda, cpg1_lim)
          : kmin(1.0f + cpg2 * lambda + cpg3 * kmin(lambda + 1.0f, 0.0f), cpg2_lim);
  return kmax(f, 0.0f);
}

// dV/dy = n . grad(n . U) with n the unit wall-normal direction (TMR's general
// geometry form). Tensor rows are velocity components, as in sst_rates.
KUTTA_KERNEL void transition_normal_derivative(int i,
    KIn(float) nx, KIn(float) ny, KIn(float) nz,
    KIn(float) ux, KIn(float) uy, KIn(float) uz,
    KIn(float) vx, KIn(float) vy, KIn(float) vz,
    KIn(float) wx, KIn(float) wy, KIn(float) wz,
    KOut(float) dvdy) {
  float a = nx[i], b = ny[i], c = nz[i];
  dvdy[i] = a * (a*ux[i] + b*uy[i] + c*uz[i])
          + b * (a*vx[i] + b*vy[i] + c*vz[i])
          + c * (a*wx[i] + b*wy[i] + c*wz[i]);
}

// Gamma sources per volume, S = source - sink*gamma (transport contract of
// sst.hpp). P - E = a*g*(1-g) - b*g*(ce2*g - 1) = (a+b)*g - (a+ce2*b)*g^2:
// the positive part stays explicit and the quadratic part becomes an implicit
// sink evaluated at the current gamma, so both coefficients are nonnegative
// and the steady source is unchanged.
KUTTA_KERNEL void transition_sources(int i,
    KIn(float) rho, KIn(float) mu, KIn(float) k, KIn(float) omega, KIn(float) distance,
    KIn(float) dvdy, KIn(float) strain, KIn(float) vorticity, KIn(float) gamma,
    KOut(float) source, KOut(float) sink, KOut(float) re_theta_c, KOut(float) onset,
    float ctu1, float ctu2, float ctu3, float cpg1, float cpg2, float cpg3,
    float cpg1_lim, float cpg2_lim, float k_floor, float omega_floor) {
  float r = rho[i], m = mu[i], d = distance[i], s = strain[i];
  float ke = kmax(k[i], k_floor), w = kmax(omega[i], omega_floor);
  // Tu_L in percent, and lambda_thetaL, Fluent 4.8.2 and Eq. 4-216/217.
  float tu = kmin(100.0f * ksqrt(2.0f * ke / 3.0f) / (w * d), 100.0f);
  float lambda = transition_lambda(dvdy[i], d, m / r);
  float fpg = transition_fpg(lambda, cpg1, cpg2, cpg3, cpg1_lim, cpg2_lim);
  float critical = ctu1 + ctu2 * kexp(-ctu3 * tu * fpg);
  // Eq. 4-215 onset functions.
  float rev = r * d * d * s / m;
  float rt = r * ke / (m * w);
  float onset1 = rev / (2.2f * critical);
  float onset2 = kmin(onset1, 2.0f);
  float rt3 = rt / 3.5f;
  float onset3 = kmax(1.0f - rt3 * rt3 * rt3, 0.0f);
  float f_onset = kmax(onset2 - onset3, 0.0f);
  float rt2 = rt / 2.0f;
  float f_turb = kexp(-(rt2 * rt2) * (rt2 * rt2));
  // P = F_length*rho*S*g*(1-g)*F_onset, F_length = 100;
  // E = c_a2*rho*Omega*g*F_turb*(c_e2*g - 1), c_a2 = .06, c_e2 = 50.
  float a = 100.0f * r * s * f_onset;
  float b = 0.06f * r * vorticity[i] * f_turb;
  float g = gamma[i];
  source[i] = (a + b) * g;
  sink[i] = (a + 50.0f * b) * g;
  re_theta_c[i] = critical;
  onset[i] = f_onset;
}

// Gamma diffusivity mu + mu_t/sigma_gamma on faces; zero flux at walls is
// structural, as for SST's k (Fluent: zero normal flux, Eq. 4-204 text).
KUTTA_KERNEL void transition_face_diffusion(int i,
    KIn(int) owner, KIn(int) neighbour, KIn(float) weight, KIn(int) wall,
    KIn(float) mu, KIn(float) mut, int n_internal, KOut(float) face, float inv_sigma) {
  int c = owner[i];
  float gc = mu[c] + mut[c] * inv_sigma;
  if (i < n_internal) {
    int n = neighbour[i]; float a = weight[i];
    face[i] = a * gc + (1.0f - a) * (mu[n] + mut[n] * inv_sigma);
  } else {
    face[i] = wall[i - n_internal] ? 0.0f : gc;
  }
}

// Inflow faces (velocity inlets, and outlets or pressure inlets with inflow)
// fix gamma at the zone's intermittency; all others extrapolate the cell value
// with zero gradient. Kinds follow fv_turbulence_boundary.
KUTTA_KERNEL void transition_boundary(int i, KIn(int) owner, KIn(int) kind,
    KIn(float) inflow, KIn(float) mass, KIn(float) gamma,
    KOut(float) boundary, KOut(int) transport_kind, int n_internal) {
  int f = i + n_internal, type = kind[i];
  int fixed = type == 1 || ((type == 2 || type == 4) && mass[f] < 0.0f);
  boundary[i] = fixed ? inflow[i] : gamma[owner[f]];
  transport_kind[i] = fixed ? 0 : 1;
}

KUTTA_KERNEL void transition_clamp(int i, KOut(float) gamma, float lo, float hi) {
  gamma[i] = kmin(kmax(gamma[i], lo), hi);
}

// SST side of the coupling. F1 = max(F1, F3), F3 = exp(-(R_y/120)^8),
// R_y = rho*d*sqrt(k)/mu, protects the laminar boundary layer's k-omega
// branch. The cell diffusivities that depend on F1 are recomputed with the
// same expressions as sst_viscosity; mu_t depends on F2 only.
KUTTA_KERNEL void sst_transition_f1(int i,
    KIn(float) rho, KIn(float) mu, KIn(float) k, KIn(float) distance, KIn(float) mut,
    KOut(float) f1, KOut(float) gamma_k, KOut(float) gamma_w, float k_floor) {
  float ry = rho[i] * distance[i] * ksqrt(kmax(k[i], k_floor)) / mu[i] / 120.0f;
  float ry2 = ry * ry, ry4 = ry2 * ry2;
  float f3 = kexp(-(ry4 * ry4));
  float blend = kmax(f1[i], f3);
  f1[i] = blend;
  gamma_k[i] = mu[i] + mut[i] * sst_inv_sigma_k(blend);
  gamma_w[i] = mu[i] + mut[i] * sst_inv_sigma_w(blend);
}

// Fluent 4-206/207: production gamma*G~_k (G~_k after SST's limiter and wall
// treatment), destruction min(max(gamma, .1), 1)*Y_k. P_k^lim adds
// 5*C_k*max(gamma-.2, 0)*(1-gamma)*F_on^lim*max(3*C_sep*mu - mu_t, 0)*S*Omega,
// F_on^lim = min(max(Re_V/(2.2*1100) - 1, 0), 3). Omega production is not
// modified (Fluent 4.7.2.2).
KUTTA_KERNEL void sst_transition_k(int i,
    KIn(float) gamma, KIn(float) rho, KIn(float) mu, KIn(float) distance,
    KIn(float) strain, KIn(float) vorticity, KIn(float) mut,
    KOut(float) k_source, KOut(float) k_sink, KOut(float) production_limit,
    float ck, float csep) {
  float g = gamma[i], d = distance[i], s = strain[i];
  float rev = rho[i] * d * d * s / mu[i];
  float on = kmin(kmax(rev / (2.2f * 1100.0f) - 1.0f, 0.0f), 3.0f);
  float extra = 5.0f * ck * kmax(g - 0.2f, 0.0f) * (1.0f - g) * on
              * kmax(3.0f * csep * mu[i] - mut[i], 0.0f) * s * vorticity[i];
  production_limit[i] = extra;
  k_source[i] = g * k_source[i] + extra;
  k_sink[i] = kmin(kmax(g, 0.1f), 1.0f) * k_sink[i];
}
