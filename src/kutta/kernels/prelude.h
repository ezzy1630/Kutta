// Kernel prelude: the restricted C++ dialect every *.kern.h is written in. The same
// kernel body compiles as plain C++ (CPU reference), CUDA C++ and Metal Shading
// Language. tools/kernelgen.py generates each backend's entry points.
//
// Kernel rules (enforced by review; kernelgen rejects bad signatures):
// - Signature: `KUTTA_KERNEL void name(int i, <params>)`. `i` is the global index;
//   bounds checking against the launch size happens in generated wrappers.
// - Params: `KOut(T) p` (read-write buffer), `KIn(T) p` (read-only buffer), or a
//   scalar `float` / `int` / `uint`. T is float, int or uint.
// - No double, no std::, no recursion, no dynamic memory, no pointer arithmetic
//   beyond indexing, no float3 in buffers (vectors are stored as separate x/y/z
//   arrays). Use the k* math wrappers below, never bare sqrt/fabs/exp.
// - Gather, don't scatter: each thread writes only its own index. No atomics, so
//   every backend is deterministic run to run.
#ifndef KUTTA_KERNEL_PRELUDE_H
#define KUTTA_KERNEL_PRELUDE_H

#if defined(KUTTA_BACKEND_METAL)
#include <metal_stdlib>
using namespace metal;
#define KUTTA_FN inline
#define KOut(T) device T*
#define KIn(T) const device T*
#define KUTTA_SQRT metal::sqrt
#define KUTTA_FMA metal::fma
#define KUTTA_ABS metal::fabs
#define KUTTA_MAX metal::fmax
#define KUTTA_MIN metal::fmin
#define KUTTA_EXP metal::exp
#define KUTTA_LOG metal::log
#define KUTTA_TANH metal::tanh
#define KUTTA_POW metal::pow
#elif defined(__CUDACC__)
#define KUTTA_FN __device__ __forceinline__
#define KOut(T) T* __restrict__
#define KIn(T) const T* __restrict__
typedef unsigned int uint;
#define KUTTA_SQRT sqrtf
#define KUTTA_FMA fmaf
#define KUTTA_ABS fabsf
#define KUTTA_MAX fmaxf
#define KUTTA_MIN fminf
#define KUTTA_EXP expf
#define KUTTA_LOG logf
#define KUTTA_TANH tanhf
#define KUTTA_POW powf
#else
#include <cmath>
#define KUTTA_FN inline
#define KOut(T) T* __restrict
#define KIn(T) const T* __restrict
typedef unsigned int uint;
#define KUTTA_SQRT std::sqrt
#define KUTTA_FMA std::fma
#define KUTTA_ABS std::fabs
#define KUTTA_MAX std::fmax
#define KUTTA_MIN std::fmin
#define KUTTA_EXP std::exp
#define KUTTA_LOG std::log
#define KUTTA_TANH std::tanh
#define KUTTA_POW std::pow
#endif

#define KUTTA_KERNEL KUTTA_FN

KUTTA_FN float kfma(float a, float b, float c) { return KUTTA_FMA(a, b, c); }
KUTTA_FN float ksqrt(float x) { return KUTTA_SQRT(x); }
KUTTA_FN float kabs(float x) { return KUTTA_ABS(x); }
KUTTA_FN float kmax(float a, float b) { return KUTTA_MAX(a, b); }
KUTTA_FN float kmin(float a, float b) { return KUTTA_MIN(a, b); }
KUTTA_FN float kexp(float x) { return KUTTA_EXP(x); }
KUTTA_FN float klog(float x) { return KUTTA_LOG(x); }
KUTTA_FN float ktanh(float x) { return KUTTA_TANH(x); }
KUTTA_FN float kpow(float x, float y) { return KUTTA_POW(x, y); }

struct kvec3 {
  float x, y, z;
};
KUTTA_FN kvec3 kv3(float x, float y, float z) {
  kvec3 v;
  v.x = x;
  v.y = y;
  v.z = z;
  return v;
}
KUTTA_FN kvec3 operator+(kvec3 a, kvec3 b) { return kv3(a.x + b.x, a.y + b.y, a.z + b.z); }
KUTTA_FN kvec3 operator-(kvec3 a, kvec3 b) { return kv3(a.x - b.x, a.y - b.y, a.z - b.z); }
KUTTA_FN kvec3 operator*(float s, kvec3 a) { return kv3(s * a.x, s * a.y, s * a.z); }
KUTTA_FN float kdot(kvec3 a, kvec3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
KUTTA_FN float knorm(kvec3 a) { return ksqrt(kdot(a, a)); }

#endif
