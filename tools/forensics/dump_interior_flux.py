"""Write Fluent's saved interior face mass flux as "c0 c1 SV_FLUX" rows (0-based cells).

SV_FLUX is positive from c0 to c1. bench_v6_operators reads this file to compare
Kutta's mass predictor with Fluent's saved flux at the saved state.
Usage: python -I dump_interior_flux.py case_dir stem out.txt   (needs numpy, h5py; see cff.py)
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np

import cff


def main():
    case_dir, stem, out_path = sys.argv[1:4]
    mesh = cff.read_mesh(os.path.join(case_dir, stem + ".cas.h5"))
    faces = cff.read_fields(os.path.join(case_dir, stem + ".dat.h5"), "faces")
    flux = cff.cell_field(faces, "SV_FLUX", mesh["c0"].size)
    interior = np.nonzero(mesh["c1"] >= 0)[0]
    if np.isnan(flux[interior]).any():
        raise SystemExit("SV_FLUX missing on interior faces")
    with open(out_path, "w") as out:
        out.write(f"{interior.size} {mesh['n_cells']}\n")
        for f in interior:
            out.write(f"{mesh['c0'][f]} {mesh['c1'][f]} {flux[f]:.17g}\n")
    print(f"interior faces {interior.size}, cells {mesh['n_cells']}")


if __name__ == "__main__":
    main()
