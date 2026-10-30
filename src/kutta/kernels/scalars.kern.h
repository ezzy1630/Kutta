#include "kutta/kernels/prelude.h"

// Device scalar wire format: normalized (hi, lo) floats, about 48 significant
// bits. Division retains the low words until the final float coefficient.
struct kpair { float hi, lo; };
KUTTA_FN kpair kp(float hi, float lo) { kpair p; p.hi = hi; p.lo = lo; return p; }
KUTTA_FN kpair kp_load(KIn(float) p, int offset) { return kp(p[offset], p[offset + 1]); }
KUTTA_FN kpair kp_add(kpair a, kpair b) {
  float s = a.hi + b.hi;
  if (s - s != 0.0f) return kp(s, 0.0f);
  float v = s - a.hi;
  float e = (a.hi - (s - v)) + (b.hi - v);
  e = e + a.lo + b.lo;
  float h = s + e;
  return kp(h, e - (h - s));
}
KUTTA_FN kpair kp_neg(kpair a) { return kp(-a.hi, -a.lo); }
KUTTA_FN kpair kp_mul(kpair a, kpair b) {
  float h = a.hi * b.hi;
  float l = kfma(a.hi, b.hi, -h) + a.hi * b.lo + a.lo * b.hi;
  return kp_add(kp(h, 0.0f), kp(l, 0.0f));
}
KUTTA_FN kpair kp_div(kpair a, kpair b) {
  float q = a.hi / b.hi;
  kpair r = kp_add(a, kp_neg(kp_mul(kp(q, 0.0f), b)));
  return kp_add(kp(q, 0.0f), kp((r.hi + r.lo) / b.hi, 0.0f));
}
KUTTA_FN void kp_store(KOut(float) out, int offset, kpair p) {
  out[offset] = p.hi; out[offset + 1] = p.lo;
}
KUTTA_FN float kp_float(kpair a) { return a.hi + a.lo; }
KUTTA_FN bool kp_positive(kpair a) { return a.hi > 0.0f || (a.hi == 0.0f && a.lo > 0.0f); }
KUTTA_FN bool kp_finite(kpair a) { return a.hi - a.hi == 0.0f && a.lo - a.lo == 0.0f; }
KUTTA_FN bool kp_zero(kpair a) { return a.hi == 0.0f && a.lo == 0.0f; }
KUTTA_FN bool kp_le(kpair a, kpair b) {
  return a.hi < b.hi || (a.hi == b.hi && a.lo <= b.lo);
}
KUTTA_FN bool kp_coefficient(kpair a) {
  float value = kp_float(a);
  // A nonzero pair that rounds to zero is an unrepresentable coefficient.
  return kp_finite(a) && value - value == 0.0f && (value != 0.0f || kp_zero(a));
}

// Scalar slots: rho=0, pAp/rhat.v/tt=2, next-rho/ts=4, z.Ap=6,
// alpha=8, beta=10, omega=12, previous-rho=14, s.s=16, r.r=18,
// tolerance^2=20, divergence-limit^2=22.
// Status: phase (0 running, 1 reliable residual requested, 2 breakdown),
// attempted updates, and the BiCGStab s-only exit. Only scalar kernels write it.
KUTTA_KERNEL void solver_reset(int i, KOut(float) scalars, KOut(int) status, int iterations) {
  if (i != 0) return;
  status[0] = 0; status[1] = iterations; status[2] = 0;
  scalars[8] = 1.0f; scalars[9] = 0.0f;
  scalars[12] = 1.0f; scalars[13] = 0.0f;
  scalars[14] = 1.0f; scalars[15] = 0.0f;
}
KUTTA_KERNEL void solver_pcg_alpha(int i, KOut(float) scalars, KOut(int) status) {
  if (i != 0) return;
  if (status[0] != 0) return;
  kpair rho = kp_load(scalars, 0), pq = kp_load(scalars, 2);
  kpair alpha = kp_div(rho, pq);
  if (!kp_positive(rho) || !kp_positive(pq) || !kp_coefficient(alpha) || kp_zero(alpha))
    status[0] = 2;
  else kp_store(scalars, 8, alpha);
}
KUTTA_KERNEL void solver_pcg_beta(int i, KOut(float) scalars, KOut(int) status, int flexible) {
  if (i != 0) return;
  if (status[0] != 0) return;
  kpair next = kp_load(scalars, 4);
  kpair numer = flexible != 0 ? kp_neg(kp_mul(kp_load(scalars, 8), kp_load(scalars, 6))) : next;
  kpair beta = kp_div(numer, kp_load(scalars, 0));
  if (!kp_positive(next) || !kp_coefficient(beta) || (kp_zero(beta) && !kp_zero(numer))) status[0] = 2;
  else { kp_store(scalars, 10, beta); kp_store(scalars, 0, next); }
}
KUTTA_KERNEL void solver_bicg_beta(int i, KOut(float) scalars, KOut(int) status) {
  if (i != 0) return;
  if (status[0] != 0) return;
  kpair rho = kp_load(scalars, 0), omega = kp_load(scalars, 12);
  kpair beta = kp_mul(kp_div(rho, kp_load(scalars, 14)), kp_div(kp_load(scalars, 8), omega));
  if (kp_zero(rho) || !kp_finite(rho) || kp_zero(omega) || !kp_coefficient(omega) || !kp_coefficient(beta) || kp_zero(beta))
    status[0] = 2;
  else { kp_store(scalars, 10, beta); status[2] = 0; }
}
KUTTA_KERNEL void solver_bicg_alpha(int i, KOut(float) scalars, KOut(int) status) {
  if (i != 0) return;
  if (status[0] != 0) return;
  kpair rv = kp_load(scalars, 2);
  kpair alpha = kp_div(kp_load(scalars, 0), rv);
  if (kp_zero(rv) || !kp_finite(rv) || !kp_coefficient(alpha) || kp_zero(alpha)) status[0] = 2;
  else kp_store(scalars, 8, alpha);
}
KUTTA_KERNEL void solver_bicg_s(int i, KIn(float) scalars, KOut(int) status) {
  if (i != 0) return;
  if (status[0] != 0) return;
  ++status[1];
  status[2] = kp_le(kp_load(scalars, 16), kp_load(scalars, 20)) ? 1 : 0;
}
KUTTA_KERNEL void solver_bicg_omega(int i, KOut(float) scalars, KOut(int) status) {
  if (i != 0) return;
  if (status[0] != 0) return;
  if (status[2] != 0) { kp_store(scalars, 12, kp(0.0f, 0.0f)); return; }
  kpair tt = kp_load(scalars, 2), ts = kp_load(scalars, 4);
  kpair omega = kp_div(ts, tt);
  if (!kp_positive(tt) || !kp_finite(tt) || kp_zero(omega) || !kp_coefficient(omega)) status[0] = 2;
  else kp_store(scalars, 12, omega);
}
KUTTA_KERNEL void solver_check(int i, KOut(float) scalars, KOut(int) status, int pcg) {
  if (i != 0) return;
  if (status[0] != 0) return;
  if (pcg != 0) ++status[1];
  kpair rr = kp_load(scalars, 18);
  if (!kp_finite(rr) || kp_le(rr, kp_load(scalars, 20)) ||
      !kp_le(rr, kp_load(scalars, 22)) || status[1] % 50 == 0) status[0] = 1;
  if (pcg == 0) kp_store(scalars, 14, kp_load(scalars, 0));
}
KUTTA_KERNEL void solver_step(int i, KOut(float) x, KOut(float) r, KIn(float) p,
                              KIn(float) q, KIn(float) scalars, KIn(int) status) {
  if (status[0] != 0) return;
  float alpha = kp_float(kp_load(scalars, 8));
  x[i] += alpha * p[i]; r[i] -= alpha * q[i];
}
KUTTA_KERNEL void solver_direction(int i, KIn(float) z, KOut(float) p,
                                   KIn(float) scalars, KIn(int) status) {
  if (status[0] == 0) p[i] = z[i] + kp_float(kp_load(scalars, 10)) * p[i];
}
KUTTA_KERNEL void solver_bicg_direction(int i, KIn(float) r, KIn(float) v,
                                        KOut(float) p, KIn(float) scalars, KIn(int) status) {
  if (status[0] != 0) return;
  float beta = kp_float(kp_load(scalars, 10)), omega = kp_float(kp_load(scalars, 12));
  p[i] = r[i] + beta * (p[i] - omega * v[i]);
}
KUTTA_KERNEL void solver_bicg_difference(int i, KIn(float) r, KIn(float) q,
                                         KOut(float) s, KIn(float) scalars, KIn(int) status) {
  if (status[0] == 0) s[i] = r[i] - kp_float(kp_load(scalars, 8)) * q[i];
}
KUTTA_KERNEL void solver_bicg_step(int i, KOut(float) x, KOut(float) r,
                                   KIn(float) z, KIn(float) zs, KIn(float) s,
                                   KIn(float) t, KIn(float) scalars, KIn(int) status) {
  if (status[0] != 0) return;
  float alpha = kp_float(kp_load(scalars, 8)), omega = kp_float(kp_load(scalars, 12));
  if (status[2] != 0) { x[i] += alpha * z[i]; r[i] = s[i]; }
  else { x[i] += alpha * z[i] + omega * zs[i]; r[i] = s[i] - omega * t[i]; }
}
KUTTA_KERNEL void solver_project(int i, KOut(float) x, KIn(float) sum, float count_hi, float count_lo) {
  x[i] -= kp_float(kp_div(kp_load(sum, 0), kp(count_hi, count_lo)));
}

// Optional K-cycle's two inner corrections also use device scalars. Phase 2
// selects its original first-correction fallback; phase 1 retains the first
// accepted correction (adaptive exit or a later breakdown).
KUTTA_KERNEL void amg_inner_alpha(int i, KOut(float) scalars, KOut(int) status,
                                  KIn(int) parent, int first) {
  if (i != 0) return;
  if (first != 0) status[0] = parent[0] == 0 ? 0 : 2;
  if (status[0] != 0) return;
  kpair denom = kp_load(scalars, 0), alpha = kp_div(kp_load(scalars, 2), denom);
  if (!kp_positive(denom) || !kp_coefficient(alpha) ||
      (kp_zero(alpha) && !kp_zero(kp_load(scalars, 2)))) status[0] = first != 0 ? 2 : 1;
  else kp_store(scalars, 8, alpha);
}
KUTTA_KERNEL void amg_inner_first(int i, KIn(float) b, KIn(float) p, KIn(float) q,
                                  KOut(float) x, KOut(float) t, KIn(float) scalars, KIn(int) status) {
  float alpha = kp_float(kp_load(scalars, 8));
  x[i] = status[0] == 0 ? alpha * p[i] : p[i];
  t[i] = status[0] == 0 ? b[i] - alpha * q[i] : 0.0f;
}
KUTTA_KERNEL void amg_inner_check(int i, KIn(float) scalars, KOut(int) status) {
  if (i != 0 || status[0] != 0) return;
  if (kp_le(kp_load(scalars, 4), kp_mul(kp(0.0625f, 0.0f), kp_load(scalars, 6)))) status[0] = 1;
}
KUTTA_KERNEL void amg_inner_beta(int i, KOut(float) scalars, KOut(int) status) {
  if (i != 0 || status[0] != 0) return;
  kpair beta = kp_div(kp_neg(kp_load(scalars, 4)), kp_load(scalars, 0));
  if (!kp_coefficient(beta) || (kp_zero(beta) && !kp_zero(kp_load(scalars, 4)))) status[0] = 1;
  else kp_store(scalars, 10, beta);
}
KUTTA_KERNEL void amg_inner_step(int i, KOut(float) x, KIn(float) p,
                                 KIn(float) scalars, KIn(int) status) {
  if (status[0] == 0) x[i] += kp_float(kp_load(scalars, 8)) * p[i];
}
