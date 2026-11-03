#!/usr/bin/env python3
"""Check a convergence benchmark's emitted controls before accepting its data.

The expected JSON object specifies numeric CASE fields such as
{"compensated_momentum": 0, "paired_state_slack": 1}. This validates the
executed binary rather than relying on its filename or build directory.
It does not compare numerical assertions or decide convergence acceptance.
"""
import argparse
import json
import math
from pathlib import Path


def check_controls(lines, expected):
    if not isinstance(expected, dict) or not expected:
        raise ValueError("expected controls must be a nonempty JSON object")
    for key, value in expected.items():
        if not isinstance(key, str) or type(value) not in (int, float):
            raise ValueError("expected control names must be strings and values numeric")
        if not math.isfinite(value):
            raise ValueError(f"expected control {key} must be finite")

    cases = [line for line in lines if line.startswith("CASE ")]
    if len(cases) != 1:
        raise ValueError(f"expected exactly one emitted CASE, found {len(cases)}")
    fields = {}
    for word in cases[0].split()[1:]:
        key, separator, value = word.partition("=")
        if not separator or key in fields:
            raise ValueError(f"invalid or duplicate CASE field: {word}")
        fields[key] = value
    for key, value in expected.items():
        if key not in fields:
            raise ValueError(f"emitted CASE is missing control {key}")
        actual = float(fields[key])
        if not math.isfinite(actual) or actual != value:
            raise ValueError(f"control {key}: expected {value}, binary emitted {fields[key]}")
    return cases[0]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--expected", type=Path, required=True)
    parser.add_argument("logs", type=Path, nargs="+")
    args = parser.parse_args()
    try:
        expected = json.loads(args.expected.read_text())
        for log in args.logs:
            check_controls(log.read_text().splitlines(), expected)
    except (OSError, ValueError) as error:
        parser.error(str(error))
    print(f"Verified emitted controls in {len(args.logs)} benchmark logs")


if __name__ == "__main__":
    main()
