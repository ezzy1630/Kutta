#!/usr/bin/env python3
"""Coupled vs SIMPLEC fields from two bench_fsae_aero --cells runs on one mesh.

Both runs renumber the same mesh identically, so cells.csv and walls.csv rows
correspond. Prints volume-weighted relative L2 differences (norm of the
difference over the norm of the SIMPLEC field) for velocity, pressure, k and
omega, and the largest wall Cp difference on the force walls.

  python3 validation/results/coupled_external/compare_fields.py <simplec run> <coupled run> [wall]
"""
import csv
import math
import pathlib
import sys


def rows(path):
    with open(path) as f:
        return list(csv.DictReader(f))


def main(simplec, coupled, wall="wing"):
    a, b = rows(pathlib.Path(simplec) / "cells.csv"), rows(pathlib.Path(coupled) / "cells.csv")
    assert len(a) == len(b) and all(x["x"] == y["x"] and x["y"] == y["y"] for x, y in zip(a, b)), "cell order differs"
    out = {}
    for name, keys in (("velocity", ("u", "v", "w")), ("pressure", ("p",)), ("k", ("k",)), ("omega", ("omega",))):
        num = den = 0.0
        for x, y in zip(a, b):
            vol = float(x["volume"])
            num += vol * sum((float(y[k]) - float(x[k])) ** 2 for k in keys)
            den += vol * sum(float(x[k]) ** 2 for k in keys)
        out[name] = math.sqrt(num / den)
    wa, wb = rows(pathlib.Path(simplec) / "walls.csv"), rows(pathlib.Path(coupled) / "walls.csv")
    cp = max(abs(float(x["cp"]) - float(y["cp"])) for x, y in zip(wa, wb) if x["patch"] == wall)
    print(" ".join(f"{k}={v:.2e}" for k, v in out.items()) + f" max_wall_cp_diff={cp:.2e}")


if __name__ == "__main__":
    main(*sys.argv[1:])
