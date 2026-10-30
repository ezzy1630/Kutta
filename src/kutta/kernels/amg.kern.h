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
// Neumaier summation retains weak couplings when large diagonal/off-diagonal
// contributions cancel; ordinary Kahan can lose the weak term entirely.
KUTTA_KERNEL void amg_galerkin(int i, KIn(int) offsets, KIn(int) entries,
                              KIn(float) fine, KOut(float) coarse) {
  float sum = 0.0f;
  float error = 0.0f;
  for (int j = offsets[i]; j < offsets[i + 1]; ++j) {
    float term = fine[entries[j]];
    float next = sum + term;
    error += kabs(sum) >= kabs(term) ? (sum - next) + term : (term - next) + sum;
    sum = next;
  }
  coarse[i] = sum + error;
}

// Memory-saving alternative, timed by bench_galerkin. One thread owns a
// coarse row and gathers directly through its member rows; no fine-nnz map.
KUTTA_KERNEL void amg_galerkin_rows(int i, KIn(int) coarse_rows, KIn(int) coarse_cols,
                                   KIn(int) offsets, KIn(int) members, KIn(int) rows,
                                   KIn(int) cols, KIn(int) aggregate, KIn(float) fine,
                                   KOut(float) coarse) {
  for (int entry = coarse_rows[i]; entry < coarse_rows[i + 1]; ++entry) {
    float sum = 0.0f;
    float error = 0.0f;
    for (int m = offsets[i]; m < offsets[i + 1]; ++m) {
      int row = members[m];
      for (int j = rows[row]; j < rows[row + 1]; ++j) {
        if (aggregate[cols[j]] == coarse_cols[entry]) {
          float term = fine[j];
          float next = sum + term;
          error += kabs(sum) >= kabs(term) ? (sum - next) + term : (term - next) + sum;
          sum = next;
        }
      }
    }
    coarse[entry] = sum + error;
  }
}

KUTTA_KERNEL void amg_restrict_guarded(int i, KIn(int) offsets, KIn(int) members,
                                      KIn(float) fine, KOut(float) coarse, KIn(int) status) {
  if (status[0] == 0) amg_restrict(i, offsets, members, fine, coarse);
}
KUTTA_KERNEL void amg_prolong_guarded(int i, KIn(int) aggregate, KIn(float) coarse,
                                     KOut(float) fine, KIn(int) status) {
  if (status[0] == 0) amg_prolong(i, aggregate, coarse, fine);
}
