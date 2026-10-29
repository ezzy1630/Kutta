"""Print the turbulence / numerics rpvars of a .cas.h5 relevant to SST replication.

Text-only scan of the Scheme strings; nothing is evaluated.  Duplicate rpvar
entries (e.g. `(k/scheme 1)` and a typed 5-element variant) are all shown.
"""
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import h5py

import cff

PATTERNS = [
    r"^kw-", r"^sst-", r"^rp-(kw|ke|sa|trans|turb|lam|visc|seg|3d|double|unsteady)", r"^turb-viscosity/",
    r"^(k|omega|mom|pressure|flow)/scheme$", r"^(k|omega|mom)/(relax|pseudo-relax)$", r"^recon/(cell-lsq|cell-lsf|lscb-wf|limit|lim/type)",
    r"^pbs/flux", r"^wall-distance", r"^curvature-correction", r"^turb-wf-", r"^sa-kappa$", r"^ke-elog-udf$", r"^ke-yplus-udf$",
    r"^pseudo-time-method/coupled-pbns/dt-method$", r"^pseudo-auto-time-step", r"^number-of-iterations$",
]


def main():
    case_dir, stem, _ = cff.args()
    with h5py.File(os.path.join(case_dir, stem + ".cas.h5"), "r") as f:
        print("Origin:", f["settings/Origin"][0].decode(), "| Solver:", f["settings/Solver"][0].decode())
        text = f["settings/Rampant Variables"][0].decode()
    flat = text.replace("\n", " ")
    seen = []
    for mm in re.finditer(r"\(([A-Za-z0-9/_\-?*]+) ((?:[^()]|\([^()]*\)){0,200})\)", flat):
        name, val = mm.group(1), mm.group(2).strip()
        if any(re.search(p, name) for p in PATTERNS) and (name, val) not in seen:
            seen.append((name, val))
    for name, val in sorted(seen):
        print(f"{name} = {val}")


if __name__ == "__main__":
    main()
