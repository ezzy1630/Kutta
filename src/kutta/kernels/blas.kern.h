// Level-1 vector kernels used by the Krylov solvers.
#include "kutta/kernels/prelude.h"

// y = y + a * x
KUTTA_KERNEL void axpy(int i, KOut(float) y, KIn(float) x, float a) { y[i] = y[i] + a * x[i]; }

// y = x + b * y
KUTTA_KERNEL void xpby(int i, KOut(float) y, KIn(float) x, float b) { y[i] = x[i] + b * y[i]; }

// z = a * x + b * y
KUTTA_KERNEL void axpbyz(int i, KOut(float) z, KIn(float) x, KIn(float) y, float a, float b) {
  z[i] = a * x[i] + b * y[i];
}

// y = a * x
KUTTA_KERNEL void scale_copy(int i, KOut(float) y, KIn(float) x, float a) { y[i] = a * x[i]; }

// z = x * y (elementwise), e.g. Jacobi preconditioner application
KUTTA_KERNEL void mul(int i, KOut(float) z, KIn(float) x, KIn(float) y) { z[i] = x[i] * y[i]; }
