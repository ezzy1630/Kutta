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
#define KUTTA_NEXTAFTER metal::nextafter
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
#define KUTTA_NEXTAFTER nextafterf
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
#define KUTTA_NEXTAFTER std::nextafter
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

// TwoSum exposes the exact rounding error of a finite float addition. Advance
// only when round-to-nearest went down; unlike an epsilon shift this retains
// exactly representable sums, including zero. No FP64 is needed on Metal.
KUTTA_FN float kadd_up(float a, float b) {
  float sum = a + b;
  float part = sum - a;
  float error = (a - (sum - part)) + (b - part);
  return error > 0.0f ? KUTTA_NEXTAFTER(sum, INFINITY) : sum;
}

// A tight upper bound for a sum: TwoSum retains the high-word error and only
// the low-word accumulation rounds upward. The final ceiling usually costs
// at most one ulp, rather than one ulp for every face of a polyhedron.
struct kupper_sum { float sum, error; };
KUTTA_FN kupper_sum kupper_add(kupper_sum a, float b) {
  float sum = a.sum + b;
  float part = sum - a.sum;
  float error = (a.sum - (sum - part)) + (b - part);
  kupper_sum result = {sum, kadd_up(a.error, error)};
  return result;
}
KUTTA_FN float kupper_value(kupper_sum a) { return kadd_up(a.sum, a.error); }

// A value carried as the unevaluated float pair hi + lo, about 48 bits.
struct kpair { float hi, lo; };
// Knuth's TwoSum: a + b == s.hi + s.lo exactly, for any finite a and b.
KUTTA_FN kpair ktwo_sum(float a, float b) {
  float sum = a + b;
  float part = sum - a;
  kpair s = {sum, (a - (sum - part)) + (b - part)};
  return s;
}
// FMA retains the product error; TwoSum carries the reduction error.
KUTTA_FN kpair kpair_product_add(kpair a,float x,float y) {
  float product=x*y;
  kpair sum=ktwo_sum(a.hi,product);
  return ktwo_sum(sum.hi,(a.lo+sum.lo)+kfma(x,y,-product));
}

KUTTA_FN kpair kpair_trace(float p,float low,float gx,float gy,float gz,
    float lx,float ly,float lz,float x,float y,float z) {
  kpair value=ktwo_sum(p,low);
  value=kpair_product_add(value,gx,x);value=kpair_product_add(value,lx,x);
  value=kpair_product_add(value,gy,y);value=kpair_product_add(value,ly,y);
  value=kpair_product_add(value,gz,z);value=kpair_product_add(value,lz,z);
  return value;
}
// Difference without first rounding either carried value to one float.
KUTTA_FN float kpair_delta(kpair a,kpair b) {
  kpair d=ktwo_sum(a.hi,-b.hi);
  return d.hi+(d.lo+(a.lo-b.lo));
}
// Adds d to the pair. hi moves only once lo exceeds slack, one ulp of the
// scale of what hi feeds (a cell value alone is not that scale: a tiny
// flux or velocity component still re-rounds the sums it enters). Near
// convergence every float quantity derived from hi (coefficients, fluxes,
// gradients) is then fixed, so lo converges against a fixed operator. A
// rounding-to-nearest hi keeps flipping at rounding midpoints, and each
// flip kicks the solution by about one ulp
// (validation/results/paired_state/REPORT.md).
KUTTA_FN kpair kpair_add(float hi, float lo, float d, float slack) {
  float low = lo + d;
  if (kabs(low) <= slack) {
    kpair kept = {hi, low};
    return kept;
  }
  kpair s = ktwo_sum(hi, d);
  return ktwo_sum(s.hi, s.lo + lo);
}

// Normalize overlapping words before adding a correction at the local
// pressure scale; global velocity-based slack can erase such corrections.
KUTTA_FN kpair kpair_add_relative(float hi,float lo,float d,float relative_slack) {
  kpair base=ktwo_sum(hi,lo),increment=ktwo_sum(base.lo,d);
  kpair value=ktwo_sum(base.hi,increment.hi);
  value=ktwo_sum(value.hi,value.lo+increment.lo);
  float slack=relative_slack*kmax(kabs(base.hi),kabs(value.hi));
  kpair difference=ktwo_sum(value.hi,-hi);
  float offset=difference.hi+(difference.lo+value.lo);
  if(kabs(offset)<=slack) {kpair kept={hi,offset};return kept;}
  return value;
}
// A held flux high word must retain the upwind sign of the physical flux.
KUTTA_FN kpair kpair_add_flux(float hi,float lo,float d,float slack) {
  kpair base=ktwo_sum(hi,lo),increment=ktwo_sum(base.lo,d);
  kpair value=ktwo_sum(base.hi,increment.hi);
  value=ktwo_sum(value.hi,value.lo+increment.lo);
  kpair difference=ktwo_sum(value.hi,-hi);
  float offset=difference.hi+(difference.lo+value.lo);
  int same=(hi>0.0f && value.hi>0.0f) || (hi<0.0f && value.hi<0.0f) || (hi==0.0f && value.hi==0.0f);
  if(same && kabs(offset)<=slack) {kpair kept={hi,offset};return kept;}
  return value;
}

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
