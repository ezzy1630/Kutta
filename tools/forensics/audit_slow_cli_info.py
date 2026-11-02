#!/usr/bin/env python3
"""Audit two CLI cases whose doctest INFO contains addresses and wall times.

This supplements compare_slow_assertions.py; it does not replace its raw result.
Only INFO metadata is normalized. Assertion locations, expressions, outcomes and
expanded values remain exact, including every solver value in the CLI logs.
"""
import argparse
import collections
import json
import re

from compare_slow_assertions import compare, read_assertions


CASES = {
    'Kutta cases converted from Fluent cases solve to the same bytes on every backend',
    'Native FSAE Fluent physics convert and solve to byte-identical Kutta outputs',
}
NUMBER = r'[-+]?\d+\.\d+e[-+]\d+'


def normalize_info(text, pointers):
    if re.fullmatch(r'0x[0-9a-fA-F]{9,16}', text):
        return pointers.setdefault(text, f'<pointer-{len(pointers)+1}>')
    # MeshTemp's process-specific directory, preserving the output filename.
    text = re.sub(r'/var/folders/[^ ]+/T/kutta-mesh-\d+-\d+/', '<MeshTemp>/', text)
    # A CLI table has one seconds column after its physical residual columns.
    header = re.search(r'\biter ((?:[a-z-]+ )+)sec/iter ', text)
    if header:
        residual_columns = len(header[1].split())
        rows = rf'(?<!\S)(\d+ )((?:{NUMBER} ){{{residual_columns}}})({NUMBER})(?= |$)'
        start = header.end()
        stop = text.find('Reached the iteration limit', start)
        if stop < 0:
            raise ValueError('Unrecognized CLI stopping message; do not discard the log')
        body, count = re.subn(rows, r'\1\2<seconds>', text[start:stop])
        if count != 3:
            raise ValueError(f'Expected the unchanged three-iteration CLI case, got {count} rows')
        text = text[:start] + body + text[stop:]
    text = re.sub(r'(\(\d+ iterations, )' + NUMBER + r'( s solve\))', r'\1<seconds>\2', text)
    return text


def normalize(cases):
    output = {}
    for name, assertions in cases.items():
        if name not in CASES:
            output[name] = assertions
            continue
        pointers, occurrences, mapped = {}, collections.Counter(), {}
        for key, value in assertions.items():
            filename, line, kind, expression, context, _ = key
            identity = (filename, line, kind, expression,
                        tuple(normalize_info(info, pointers) for info in context))
            occurrence = occurrences[identity]
            occurrences[identity] += 1
            mapped[(*identity, occurrence)] = value
        if len(mapped) != len(assertions):
            raise ValueError('Normalization lost an assertion')
        output[name] = mapped
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('baseline')
    parser.add_argument('branch')
    args = parser.parse_args()
    before, after = read_assertions(args.baseline), read_assertions(args.branch)
    raw = compare(before, after)
    audited = compare(normalize(before), normalize(after))
    keys = ('regressions', 'recoveries', 'changed_values', 'missing_assertions',
            'new_assertions', 'missing_cases', 'new_cases')
    result = {'scope': sorted(CASES), 'normalization': 'INFO addresses, MeshTemp paths and CLI wall-clock fields only; expanded values unchanged',
              'raw': {key: len(raw[key]) for key in keys},
              'audited': {key: len(audited[key]) for key in keys},
              'cases': audited['cases']}
    print(json.dumps(result, indent=2))
    return int(any(audited[key] for key in keys))


if __name__ == '__main__':
    raise SystemExit(main())
