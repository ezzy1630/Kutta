// Kernels that exist only so tests can check backend math accuracy.
#include "kutta/kernels/prelude.h"

// out = f(x) for f selected by fn: 0 sqrt|x|, 1 exp, 2 log(|x|+1e-3), 3 tanh, 4 (|x|+0.5)^1.7
KUTTA_KERNEL void test_math(int i, KOut(float) out, KIn(float) x, int fn) {
  const float v = x[i];
  float r;
  if (fn == 0) r = ksqrt(kabs(v));
  else if (fn == 1) r = kexp(v);
  else if (fn == 2) r = klog(kabs(v) + 1e-3f);
  else if (fn == 3) r = ktanh(v);
  else r = kpow(kabs(v) + 0.5f, 1.7f);
  out[i] = r;
}
