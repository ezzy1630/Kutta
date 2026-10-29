#!/usr/bin/env python3
"""Check Re=1000 refinement against Erturk's 601^2 centreline tables.

Generate the 65^2, 129^2 and 257^2 profile CSVs with bench_fv. Printed-value
rounding intervals are not an error tolerance for a finite-volume grid.
The doctest slow suite also reports primary-vortex Richardson estimates.
"""

import argparse
import csv
import math
from pathlib import Path


def rows(path):
    with path.open() as source:
        return list(csv.DictReader(source))


def interpolate(profile, coordinate):
    for left, right in zip(profile, profile[1:]):
        if left[0] <= coordinate <= right[0]:
            weight = (coordinate - left[0]) / (right[0] - left[0])
            return left[1] + weight * (right[1] - left[1])
    raise ValueError(f"Coordinate outside profile: {coordinate}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profiles", type=Path, default=Path("build"))
    parser.add_argument("--data", type=Path, default=Path("validation/data"))
    parser.add_argument("--backend", choices=("cpu", "metal"), action="append")
    args = parser.parse_args()
    failed = False
    for backend in args.backend or ("cpu", "metal"):
        profiles = {n: rows(args.profiles / f"fv-cavity-{backend}-Re1000-{n}x{n}.csv")
                    for n in (65, 129, 257)}
        for component, axis, direction in (("u", "y", "vertical"), ("v", "x", "horizontal")):
            reference = rows(args.data / "erturk2005_cavity" /
                             f"erturk2005_{component}_{direction}_centerline_Re1000.csv")
            errors, centers = [], []
            for n, profile in profiles.items():
                values = sorted((float(r["coordinate"]), float(r["velocity"]))
                                for r in profile if r["component"] == component)
                error = max(abs(interpolate(values, float(r[axis])) - float(r[f"{component}_Re1000"]))
                            for r in reference)
                errors.append(error)
                centers.append(interpolate(values, .5))
                print(f"CAVITY_REFINEMENT backend={backend} component={component} n={n} "
                      f"max_error={error:.12g} center={centers[-1]:.12g}")
            order = math.log(errors[1] / errors[2]) / math.log(257 / 129)
            passed = errors[2] < errors[1] and order >= 1.6
            failed |= not passed
            richardson = centers[2] + (centers[2] - centers[1]) / ((257 / 129)**2 - 1)
            print(f"CAVITY_GATE backend={backend} component={component} "
                  f"order={order:.12g} richardson_center={richardson:.12g} pass={int(passed)}")
    return int(failed)


if __name__ == "__main__":
    raise SystemExit(main())
