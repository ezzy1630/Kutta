#!/usr/bin/env python3
"""Supplement the canonical comparison for two CLI tests with volatile INFO.

Usage: python3 compare_cli_context.py BASELINE_XML BRANCH_XML > comparison.json

Run tools/forensics/compare_slow_assertions.py first and retain its result.
This audit uses that tool's assertion reader and comparator unchanged. Only
the two named CLI cases normalize standalone pointer INFO, temporary mesh
directories, the last (time) column of printed residual rows, and solve time.
All expressions, outcomes, expanded values, residuals and other INFO remain.
No numerical tolerance or assertion is changed or omitted.
"""
import importlib.util
import json
import pathlib
import re
import sys


CLI_CASES = {
    'Kutta cases converted from Fluent cases solve to the same bytes on every backend',
    'Native FSAE Fluent physics convert and solve to byte-identical Kutta outputs',
}


def normalize(data):
    for name, assertions in data.items():
        if name not in CLI_CASES:
            continue
        normalized, pointers = {}, {}

        def context(text):
            if re.fullmatch(r'0x[0-9a-fA-F]{9,16}', text):
                return pointers.setdefault(text, f'<pointer-{len(pointers) + 1}>')
            text = re.sub(r'/var/folders/[^ ]+/kutta-mesh-[0-9]+-[0-9]+',
                          '<temporary-mesh-directory>', text)
            # CLI rows contain an integer iteration and scientific-notation
            # residuals followed by seconds/iteration. Keep every residual.
            text = re.sub(
                r'(?<!\S)([0-9]+)((?: [-+]?\d+\.\d+e[-+]\d+){4,})'
                r'(?= (?:[0-9]+ |Reached the iteration limit))',
                lambda match: match[1] + ' ' + ' '.join(match[2].split()[:-1])
                + ' <seconds-per-iteration>', text)
            return re.sub(r'(\([0-9]+ iterations), [-+0-9.e]+ s solve\)',
                          r'\1, <seconds> s solve)', text)

        for key, value in assertions.items():
            identity = (*key[:4], tuple(context(info) for info in key[4]), key[5])
            if identity in normalized:
                raise ValueError(f'Normalization merged assertions in {name}')
            normalized[identity] = value
        data[name] = normalized
    return data


def main():
    tool = pathlib.Path(__file__).resolve().parents[3] / 'tools/forensics/compare_slow_assertions.py'
    spec = importlib.util.spec_from_file_location('canonical_comparison', tool)
    canonical = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(canonical)
    result = canonical.compare(normalize(canonical.read_assertions(sys.argv[1])),
                               normalize(canonical.read_assertions(sys.argv[2])))
    print(json.dumps(result, indent=2))
    return 0 if result['no_new_failures'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
