#!/usr/bin/env python3
"""Bitwise output parity of two `kutta` binaries on the memory-lane cases.

Runs `kutta solve CASE --iterations N --out DIR` for every case and backend
with each binary and compares SHA-256 of fluxes.csv, forces.csv, loads.json,
result.vtu and residuals.csv without its timing column. Prints one JSON line
per run and exits non-zero on any difference.

  python3 tools/memory_scale/parity.py --baseline BASE/kutta --candidate build/kutta \\
      --out /tmp/parity --backends cpu,metal
"""
import argparse
import hashlib
import json
import pathlib
import subprocess
import sys
import tempfile

ROOT = pathlib.Path(__file__).resolve().parents[2]
EXAMPLES = ROOT / "validation/fluent_examples/mixing_elbow"
CASES = {
    "v6_elbow": (EXAMPLES / "mixing_elbow.cas.h5", "wfgc"),
    "v6_simplec": (None, "wfgc"),  # converted by the baseline binary below
    "elbow_param": (EXAMPLES / "elbow_param.cas.h5", "compressibility,ke-options"),
    "ahmed_l1": (ROOT / "runs/mesh/ahmed35_L1.cas", "wfgc"),
    "ahmed_l1_simplec": (ROOT / "runs/mesh/ahmed35_L1_simplec.cas", "wfgc"),
    "ahmed_l1_aero_simplec": (ROOT / "tools/memory_scale/ahmed_simplec.kutta.json", ""),
}
FILES = ("fluxes.csv", "forces.csv", "loads.json", "result.vtu", "residuals.csv")


def digest(path):
    data = path.read_bytes()
    if path.name == "residuals.csv":
        # Match grep -v '^#' | cut -d, -f1-7 from the common brief.
        data = b"".join(b",".join(line.split(b",")[:7]) + b"\n" for line in data.splitlines()
                         if not line.startswith(b"#"))
    return hashlib.sha256(data).hexdigest()


def compare_files(before, after):
    """Run cmp on outputs, including the brief's literal grep/cut residual path."""
    differs = []
    with tempfile.TemporaryDirectory() as temporary:
        for name in FILES:
            paths = [before / name, after / name]
            if name == "residuals.csv":
                normalized = []
                for index, path in enumerate(paths):
                    filtered = subprocess.check_output(["grep", "-v", "^#", str(path)])
                    data = subprocess.check_output(["cut", "-d,", "-f1-7"], input=filtered)
                    target = pathlib.Path(temporary) / str(index)
                    target.write_bytes(data)
                    normalized.append(target)
                paths = normalized
            status = subprocess.run(["cmp", "-s", *(str(path) for path in paths)]).returncode
            if status > 1:
                raise RuntimeError(f"cmp failed for {name}, status {status}")
            if status:
                differs.append(name)
    return differs


def run(binary, case, approximate, backend, iterations, out, edge_storage=None):
    out.mkdir(parents=True, exist_ok=False)
    cmd = [str(binary), "solve", str(case), "--backend", backend, "--iterations", str(iterations),
           "--out", str(out)]
    if approximate:
        cmd += ["--approximate", approximate]
    if edge_storage:
        cmd += ["--coupled-edge-storage", edge_storage]
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    (out / "solve.log").write_text(proc.stdout)
    if proc.returncode:
        raise SystemExit(f"{' '.join(cmd)} failed:\n{proc.stdout[-2000:]}")
    actual = json.loads((out / "loads.json").read_text()).get("iterations")
    if actual != iterations:
        raise SystemExit(f"{case.name} on {backend} completed {actual} iterations; expected {iterations}")
    return {name: digest(out / name) for name in FILES}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--baseline", required=True)
    ap.add_argument("--candidate", required=True)
    for label in ("baseline", "candidate"):
        ap.add_argument(f"--{label}-edge-storage", choices=("cached", "recompute"))
    ap.add_argument("--out", required=True)
    ap.add_argument("--backends", default="cpu,metal")
    ap.add_argument("--cases", default=",".join(CASES))
    ap.add_argument("--iterations", type=int, default=40)
    a = ap.parse_args()
    if pathlib.Path(a.out).exists():
        ap.error("use a fresh output directory")
    same = True
    for name in a.cases.split(","):
        case, approximate = CASES[name]
        if name == "v6_simplec":
            generated = pathlib.Path(a.out) / "v6_simplec_case"
            generated.mkdir(parents=True, exist_ok=True)
            setup = generated / "setup.json"
            setup.write_text(json.dumps({"v": 1, "overrides": {"methods/pressure-velocity-coupling": "simplec"}}))
            case = generated / "simplec.kutta.json"
            subprocess.run([a.baseline, "convert", str(CASES["v6_elbow"][0]), "--setup", str(setup),
                            "--out", str(case)], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        for backend in a.backends.split(","):
            hashes = {label: run(binary, case, approximate, backend, a.iterations,
                                 pathlib.Path(a.out) / f"{name}_{backend}_{label}",
                                 getattr(a, label + "_edge_storage"))
                      for label, binary in (("baseline", a.baseline), ("candidate", a.candidate))}
            diff = [f for f in FILES if hashes["baseline"][f] != hashes["candidate"][f]]
            directories = [pathlib.Path(a.out) / f"{name}_{backend}_{label}" for label in ("baseline", "candidate")]
            diff = sorted(set(diff) | set(compare_files(*directories)))
            same &= not diff
            print(json.dumps({"case": name, "backend": backend, "iterations": a.iterations,
                              "binary_sha256": {label: hashlib.sha256(pathlib.Path(binary).read_bytes()).hexdigest()
                                                for label, binary in (("baseline", a.baseline), ("candidate", a.candidate))},
                              "case_sha256": hashlib.sha256(case.read_bytes()).hexdigest(),
                              "baseline_edge_storage": a.baseline_edge_storage or "cached",
                              "candidate_edge_storage": a.candidate_edge_storage or "cached",
                              "identical": not diff, "differs": diff, "sha256": hashes["candidate"]}), flush=True)
    sys.exit(0 if same else 1)


if __name__ == "__main__":
    main()
