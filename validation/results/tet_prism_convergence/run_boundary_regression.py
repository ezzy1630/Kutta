"""Build and run the boundary regression, rejecting stale or empty results.

The test's exit code is preserved, including during a deliberate red run.
--preserve-cli restores the existing CLI after the focused library test so
an in-progress reference run can keep its own executable.
"""
import argparse
import hashlib
import json
import pathlib
import shutil
import subprocess
import tempfile
import xml.etree.ElementTree as ET


ROOT = pathlib.Path(__file__).resolve().parents[3]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_identity():
    paths = sorted(path for path in (ROOT / 'src').rglob('*') if path.is_file())
    paths += [ROOT / 'tools/kernelgen.py', ROOT / 'tests/fv/test_native_farfield.cpp',
              ROOT / 'tests/data/meshes/native_farfield_tets.msh']
    hashes = {str(path.relative_to(ROOT)): digest(path) for path in paths}
    return hashlib.sha256(json.dumps(hashes, sort_keys=True).encode()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=pathlib.Path, help='fresh output prefix under an ignored runs folder')
    parser.add_argument('--build', type=pathlib.Path, default=ROOT / 'build')
    parser.add_argument('--case', help='additional doctest case filter, for checking the zero-case guard')
    parser.add_argument('--preserve-cli', action='store_true')
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    xml = pathlib.Path(str(args.output) + '.xml')
    if xml.exists():
        raise RuntimeError('The result XML already exists; choose a fresh output prefix')
    before = source_identity()
    cli, binary = args.build / 'kutta', args.build / 'kutta-tests'
    with tempfile.TemporaryDirectory(prefix='boundary-regression-', dir=args.output.parent) as temporary:
        backup = pathlib.Path(temporary) / 'kutta'
        if args.preserve_cli:
            shutil.copy2(cli, backup)
        try:
            cache = (args.build / 'CMakeCache.txt').read_text().splitlines()
            configured = next(line.split('=', 1)[1] for line in cache
                              if line.startswith('CMAKE_HOME_DIRECTORY:INTERNAL='))
            if pathlib.Path(configured).resolve() != ROOT:
                raise RuntimeError('The build belongs to a different checkout')
            with pathlib.Path(str(args.output) + '.build.log').open('w') as log:
                subprocess.run(['cmake', '--build', str(args.build), '--target', 'kutta-tests', '-j', '4'],
                               stdout=log, stderr=subprocess.STDOUT, check=True)
            if source_identity() != before:
                raise RuntimeError('Source changed during the build')
            command = [str(binary), '-sf=*test_native_farfield.cpp', '-r=xml', '--success']
            if args.case:
                command.append('-tc=' + args.case)
            with xml.open('w') as out:
                result = subprocess.run(command, stdout=out, check=False)
            if source_identity() != before:
                raise RuntimeError('Source changed during the test')
            root = ET.parse(xml).getroot()
            cases = [case for case in root.iter('TestCase')
                     if case.get('skipped') != 'true' and case.find('OverallResultsAsserts') is not None]
            if len(cases) != 1:
                raise RuntimeError(f'Expected exactly one executed boundary case, observed {len(cases)}')
            case = cases[0]
            if not case.get('filename', '').replace('\\', '/').endswith('/test_native_farfield.cpp'):
                raise RuntimeError('The executed case is not the boundary regression')
            contexts = [tuple(' '.join(info.itertext()).strip() for info in expression.findall('Info'))
                        for expression in case.iter('Expression')]
            combinations = {(backend, boundary) for backend in ['cpu', 'metal']
                            for boundary in ['symmetry', 'moving wall']
                            if any(f'backend: {backend}' in context and boundary in context for context in contexts)}
            if len(combinations) != 4:
                raise RuntimeError(f'Missing CPU/Metal boundary combinations: {combinations}')
            summary = case.find('OverallResultsAsserts')
            failures = int(summary.get('failures'))
            successes = int(summary.get('successes'))
            if not successes + failures or result.returncode != (1 if failures else 0):
                raise RuntimeError('Test exit disagrees with the executed assertions')
            receipt = {'source_sha256': before, 'binary_sha256': digest(binary),
                       'executed_cases': 1, 'case': case.get('name'),
                       'boundary_backend_combinations': sorted(combinations),
                       'successes': successes, 'failures': failures, 'exit': result.returncode}
            pathlib.Path(str(args.output) + '.receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
            print(json.dumps(receipt, indent=2))
            return result.returncode
        finally:
            if args.preserve_cli:
                shutil.copy2(backup, cli)
                if digest(backup) != digest(cli):
                    raise RuntimeError('The original CLI was not restored')


if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except Exception as error:
        raise SystemExit(f'Boundary regression runner rejected the result: {error}')
