"""Run the complete slow tier in separate, serial processes to bound memory.

Only doctest's first/last case range changes; no subcase, backend, solver
iteration budget, assertion, or reference is filtered. XML stays in the given
ignored output directory. The merged document retains every case assertion.
"""
import argparse
import copy
import hashlib
import json
import os
import pathlib
import re
import subprocess
import xml.etree.ElementTree as ET


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_case(xml, resources, name, returncode):
    root = ET.parse(xml).getroot()
    executed = [case for case in root.iter('TestCase')
                if case.get('skipped') != 'true' and
                case.find('OverallResultsAsserts') is not None]
    if len(executed) != 1 or executed[0].get('name') != name:
        raise RuntimeError(f'case output did not execute exactly {name}')
    case = executed[0]
    if any('crash' in error.attrib for error in case.iter('Exception')):
        raise RuntimeError(f'case {name} crashed')
    summary = case.find('OverallResultsAsserts')
    failures = int(summary.get('failures'))
    if returncode != (1 if failures else 0):
        raise RuntimeError(f'case {name} exit {returncode} disagrees with its assertions')
    match = re.search(r'(\d+)\s+peak memory footprint', resources.read_text())
    if not match:
        raise RuntimeError(f'case {name} lacks a memory readback')
    footprint = int(match[1])
    record = {'case': name, 'exit': returncode, 'failures': failures,
              'successes': int(summary.get('successes')), 'peak_footprint_bytes': footprint,
              'xml_sha256': digest(xml)}
    if footprint > 3_000_000_000:
        raise RuntimeError(f'case {name} exceeded the 3 GB envelope')
    return case, record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('binary', type=pathlib.Path)
    parser.add_argument('output', type=pathlib.Path)
    parser.add_argument('--cli', type=pathlib.Path, default=pathlib.Path('build/kutta'),
                        help='the CLI path embedded in this test build')
    parser.add_argument('--resume', action='store_true',
                        help='reuse only completed case receipts with matching binaries and XML hashes')
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    binary = args.binary.resolve()
    env = dict(os.environ, MallocSpaceEfficient='1')
    listing = subprocess.check_output(
        [str(binary), '-ts=slow', '--order-by=file', '--list-test-cases'],
        env=env, text=True)
    names = [line for line in listing.splitlines()
             if line and not line.startswith(('[doctest]', '='))]
    if len(names) != 37 or len(set(names)) != len(names):
        raise RuntimeError(f'expected the unchanged 37-case slow tier, got {len(names)}')
    cli_hash = digest(args.cli)
    manifest = {'binary_sha256': digest(binary), 'cli_sha256': cli_hash, 'cases': []}
    merged = ET.Element('doctest', {'serial_case_processes': 'true'})
    suite = ET.SubElement(merged, 'TestSuite', {'name': 'slow'})
    manifest_path = args.output / 'manifest.json'
    if args.resume:
        previous = json.loads(manifest_path.read_text())
        if any(previous.get(key) != manifest[key] for key in ['binary_sha256', 'cli_sha256']):
            raise RuntimeError('Resume binary or CLI does not match the recorded verification')
        if len(previous['cases']) > len(names):
            raise RuntimeError('Resume manifest contains too many cases')
        for number, saved in enumerate(previous['cases'], 1):
            xml = args.output / f'{number:02d}.xml'
            resources = args.output / f'{number:02d}.resources.log'
            case, record = read_case(xml, resources, names[number - 1], saved['exit'])
            if record != saved:
                raise RuntimeError(f'Resume case {number} no longer matches its receipt')
            manifest['cases'].append(record)
            suite.append(copy.deepcopy(case))
        manifest['resumed_from_cases'] = len(manifest['cases'])
        print(f'RESUME {len(manifest["cases"])}/37 immutable completed cases', flush=True)
    elif manifest_path.exists():
        raise RuntimeError('Verification output already contains receipts; use --resume')
    completed = len(manifest['cases'])
    for number, name in enumerate(names, 1):
        if number <= completed:
            continue
        if digest(args.cli) != cli_hash:
            raise RuntimeError('The test build\'s CLI changed during verification')
        print(f'START {number}/37 {name}', flush=True)
        xml = args.output / f'{number:02d}.xml'
        resources = args.output / f'{number:02d}.resources.log'
        if xml.exists() or resources.exists():
            raise RuntimeError(f'Unrecorded case {number} output exists; preserve it before resuming')
        command = ['/usr/bin/time', '-l', str(binary), '-ts=slow',
                   '--order-by=file', f'--first={number}', f'--last={number}',
                   '-r=xml', '--success']
        with xml.open('w') as out, resources.open('w') as err:
            run = subprocess.run(command, stdout=out, stderr=err, env=env)
        if digest(args.cli) != cli_hash:
            raise RuntimeError('The test build\'s CLI changed during a case')
        case, record = read_case(xml, resources, name, run.returncode)
        manifest['cases'].append(record)
        (args.output / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
        suite.append(copy.deepcopy(case))
        print(f'END {number}/37 failures={record["failures"]} footprint={record["peak_footprint_bytes"]}', flush=True)
    ET.SubElement(merged, 'OverallResultsAsserts', {
        'successes': str(sum(case['successes'] for case in manifest['cases'])),
        'failures': str(sum(case['failures'] for case in manifest['cases'])),
    })
    ET.SubElement(merged, 'OverallResultsTestCases', {
        'successes': str(sum(case['failures'] == 0 for case in manifest['cases'])),
        'failures': str(sum(case['failures'] != 0 for case in manifest['cases'])),
        'skipped': '0',
    })
    ET.indent(merged)
    ET.ElementTree(merged).write(args.output / 'slow.xml', encoding='utf-8', xml_declaration=True)
    print('COMPLETE 37/37', flush=True)


if __name__ == '__main__':
    main()
