#include "kutta/kernels/prelude.h"

// P is piecewise constant, R = P^T. Each coarse row owns its member list;
// interpolation reads the fine->aggregate map. Neither operation scatters.
KUTTA_KERNEL void amg_restrict(int i, KIn(int) offsets, KIn(int) members,
                              KIn(float) fine, KOut(float) coarse) {
  float sum = 0.0f;
  float error = 0.0f;
  for (int j = offsets[i]; j < offsets[i + 1]; ++j) {
    float term = fine[members[j]] - error;
    float next = sum + term;
    error = (next - sum) - term;
    sum = next;
  }
  coarse[i] = sum;
}

KUTTA_KERNEL void amg_prolong(int i, KIn(int) aggregate, KIn(float) coarse,
                             KOut(float) fine) {
  if (aggregate[i] >= 0) fine[i] += coarse[aggregate[i]];
}

// Precomputed inverse fine-entry -> coarse-entry map. A coarse entry owns
// its contributing fine entries, so Galerkin refresh is deterministic.
KUTTA_KERNEL void amg_galerkin(int i, KIn(int) offsets, KIn(int) entries,
                              KIn(float) fine, KOut(float) coarse) {
  float sum = 0.0f;
  float error = 0.0f;
  for (int j = offsets[i]; j < offsets[i + 1]; ++j) {
    float term = fine[entries[j]] - error;
    float next = sum + term;
    error = (next - sum) - term;
    sum = next;
  }
  coarse[i] = sum;
}
