"""Decode /results/residuals: test scaled = col0/col1 against Fluent's default 1e-3 stop criterion."""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import h5py
import numpy as np

import cff

case_dir, stem, _ = cff.args()
with h5py.File(os.path.join(case_dir, stem + ".dat.h5"), "r") as f:
    g = f["results/residuals/phase-1"]
    scaled = {}
    for eq in g:
        d = g[eq]["data"][:]
        it = g[eq]["iterations"][:]
        scaled[eq] = d[:, 0] / d[:, 1]
        print(f"{eq:11s} n={len(it)} iters {it[0]:.0f}..{it[-1]:.0f}  cols2-3 all zero: {not d[:, 2:].any()}  "
              f"col1 const after it5: {np.allclose(d[5:, 1], d[5, 1])}  col1==max(col0[:5]): {np.isclose(d[5, 1], d[:5, 0].max())}  "
              f"scaled last {scaled[eq][-1]:.3e} prev {scaled[eq][-2]:.3e}")
    worst = np.max([v for k, v in scaled.items() if k != "energy"], axis=0)
    print("first iteration where all non-energy scaled residuals < 1e-3:", int(np.argmax(worst < 1e-3)) + 1, "of", len(worst),
          "; energy < 1e-6 at end:", scaled["energy"][-1] < 1e-6)
