#!/usr/bin/env python3
"""Check published-data scatter at coordinates shared by Ghia and Erturk.

No scatter band is invented for Re=100/400, where only Ghia is supplied.
Re=1000 intervals include each table's last printed digit rounding interval.
This is a necessary common-point check, not proof of the entire profile gate.
"""

import argparse
import csv
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
    args = parser.parse_args()
    failed = False
    for backend in ("cpu", "metal"):
        profile = rows(args.profiles / f"fv-cavity-{backend}-Re1000-129x129.csv")
        for component, axis, direction in (("u", "y", "vertical"), ("v", "x", "horizontal")):
            values = sorted((float(r["coordinate"]), float(r["velocity"]))
                            for r in profile if r["component"] == component)
            ghia = rows(args.data / "ghia1982_cavity" /
                        f"ghia1982_{component}_{direction}_centerline.csv")
            erturk = rows(args.data / "erturk2005_cavity" /
                          f"erturk2005_{component}_{direction}_centerline_Re1000.csv")
            fine = {float(r[axis]): float(r[f"{component}_Re1000"]) for r in erturk}
            for point in ghia:
                coordinate = float(point[f"{axis}_exact"])
                if coordinate not in fine:
                    continue
                a, b = float(point[f"{component}_Re1000"]), fine[coordinate]
                low, high = min(a - .000005, b - .00005), max(a + .000005, b + .00005)
                actual = interpolate(values, coordinate)
                passed = low <= actual <= high
                failed |= not passed
                print(f"CAVITY_SCATTER backend={backend} component={component} "
                      f"coordinate={coordinate:g} actual={actual:.12g} "
                      f"interval=[{low:.8g},{high:.8g}] pass={int(passed)}")
    return int(failed)


if __name__ == "__main__":
    raise SystemExit(main())
