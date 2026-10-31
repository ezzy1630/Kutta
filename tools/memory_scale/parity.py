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

ROOT = pathlib.Path(__file__).resolve().parents[2]
EXAMPLES = ROOT / "validation/fluent_examples/mixing_elbow"
CASES = {
    "v6_elbow": (EXAMPLES / "mixing_elbow.cas.h5", "wfgc"),
    "elbow_param": (EXAMPLES / "elbow_param.cas.h5", "compressibility,ke-options"),
    "ahmed_l1": (ROOT / "runs/mesh/ahmed35_L1.cas", "wfgc"),
}
FILES = ("fluxes.csv", "forces.csv", "loads.json", "result.vtu", "residuals.csv")


def digest(path):
    data = path.read_bytes()
    if path.name == "residuals.csv":
        # The last column is seconds per iteration.
        data = b"\n".join(line.rsplit(b",", 1)[0] for line in data.split(b"\n"))
    return hashlib.sha256(data).hexdigest()


def run(binary, case, approximate, backend, iterations, out):
    out.mkdir(parents=True, exist_ok=True)
    cmd = [str(binary), "solve", str(case), "--backend", backend, "--iterations", str(iterations),
           "--out", str(out), "--approximate", approximate]
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    (out / "solve.log").write_text(proc.stdout)
    if proc.returncode:
        raise SystemExit(f"{' '.join(cmd)} failed:\n{proc.stdout[-2000:]}")
    return {name: digest(out / name) for name in FILES}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--baseline", required=True)
    ap.add_argument("--candidate", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--backends", default="cpu,metal")
    ap.add_argument("--cases", default=",".join(CASES))
    ap.add_argument("--iterations", type=int, default=40)
    a = ap.parse_args()
    same = True
    for name in a.cases.split(","):
        case, approximate = CASES[name]
        for backend in a.backends.split(","):
            hashes = {label: run(binary, case, approximate, backend, a.iterations,
                                 pathlib.Path(a.out) / f"{name}_{backend}_{label}")
                      for label, binary in (("baseline", a.baseline), ("candidate", a.candidate))}
            diff = [f for f in FILES if hashes["baseline"][f] != hashes["candidate"][f]]
            same &= not diff
            print(json.dumps({"case": name, "backend": backend, "iterations": a.iterations,
                              "identical": not diff, "differs": diff, "sha256": hashes["candidate"]}), flush=True)
    sys.exit(0 if same else 1)


if __name__ == "__main__":
    main()
