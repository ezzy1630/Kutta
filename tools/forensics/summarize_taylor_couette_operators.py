#!/usr/bin/env python3
"""Record every Taylor-Couette operator norm and adjacent refinement order."""

import argparse
import hashlib
import json
import math
from pathlib import Path


def summarize(path):
    groups = {}
    for line in path.read_text().splitlines():
        if not line.startswith("TC_OPERATOR "):
            continue
        fields = dict(item.split("=", 1) for item in line.split()[1:])
        level = int(fields.pop("nr"))
        rms = float(fields.pop("rms"))
        maximum = float(fields.pop("maximum"))
        key = tuple(fields.items())
        levels = groups.setdefault(key, {})
        if level in levels:
            raise ValueError(f"duplicate operator at nr={level}: {key}")
        levels[level] = (rms, maximum)
    if not groups:
        raise ValueError(f"no operator records in {path}")
    records = []
    for key, levels in groups.items():
        sizes = sorted(levels)
        errors = [levels[n][0] for n in sizes]
        orders = [
            math.log(a / b) / math.log(n / old_n)
            if 0 < a < math.inf and 0 < b < math.inf else None
            for old_n, n, a, b in zip(sizes, sizes[1:], errors, errors[1:])
        ]
        records.append(dict(key, rings=sizes, rms=errors,
                            maximum=[levels[n][1] for n in sizes], orders=orders))
    return {"input": path.name, "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "records": records}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.write_text(json.dumps(summarize(args.input), indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
