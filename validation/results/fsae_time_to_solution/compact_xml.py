#!/usr/bin/env python3
"""Publish assertion identities and value hashes without generated mesh payloads.

Compare the original --success XML with compare_slow_assertions.py first. The
originals stay local. These compact XML receipts work with the same comparator;
SHA-256 preserves value equality without republishing full fields or meshes.
"""
import argparse
import gzip
import hashlib
import importlib.util
import pathlib
import xml.etree.ElementTree as ET

ROOT = pathlib.Path(__file__).resolve().parents[3]


def digest(value):
    return 'sha256:' + hashlib.sha256(value.encode()).hexdigest()


def compact(source):
    spec = importlib.util.spec_from_file_location('slow_assertions', ROOT / 'tools/forensics/compare_slow_assertions.py')
    reader = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(reader)
    # The canonical reader checks completeness and normalizes pointer contexts.
    cases = reader.read_assertions(source)
    result = ET.Element('doctest', {'receipt': 'assertion-identities-and-value-hashes',
                                 'source_sha256': hashlib.sha256(pathlib.Path(source).read_bytes()).hexdigest()})
    for name, assertions in cases.items():
        case = ET.SubElement(result, 'TestCase', {'name': name})
        for (filename, line, kind, original, context, occurrence), value in assertions.items():
            del occurrence  # reconstructed by the canonical reader in this order
            expr = ET.SubElement(case, 'Expression', {'filename': filename, 'line': line, 'type': kind,
                                                     'success': str(value['passed']).lower()})
            original = original.replace(str(ROOT), '<checkout>')
            if len(original.encode()) > 512 or '/Users/' in original or '/home/' in original:
                original = digest(original)
            ET.SubElement(expr, 'Original').text = original
            ET.SubElement(expr, 'Expanded').text = digest(value['expanded'])
            for info in context:
                ET.SubElement(expr, 'Info').text = digest(info)
        passed = sum(x['passed'] for x in assertions.values())
        ET.SubElement(case, 'OverallResultsAsserts', {'successes': str(passed),
                                                    'failures': str(len(assertions) - passed)})
    return ET.tostring(result, encoding='utf-8', xml_declaration=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=pathlib.Path)
    parser.add_argument('output', type=pathlib.Path)
    args = parser.parse_args()
    data = compact(args.source)
    args.output.write_bytes(gzip.compress(data, mtime=0) if args.output.suffix == '.gz' else data)
