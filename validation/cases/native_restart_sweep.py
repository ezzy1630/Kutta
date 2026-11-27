"""Compare independent cold runs with forward and reverse native restarts."""
import argparse
import hashlib
import json
import os
import platform
import re
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SPEEDS = (8, 9, 10, 11, 12)


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2) + "\n")


def read_json(path):
    return json.loads(path.read_text())


def sha256(path):
    with path.open("rb") as source:
        return hashlib.file_digest(source, "sha256").hexdigest()


def torque(zones, name):
    wall = next(z for z in zones if z["name"] == name)
    return wall["viscous_moment"][2] + wall["pressure_moment"][2]


def sweep(cli, root):
    root.mkdir(parents=True)
    mesh = root / "shaft.msh"
    subprocess.run([sys.executable, ROOT / "tools/fsae_aero/rotating_shaft.py", mesh,
                    "--radial", "16", "--angular", "72"], check=True)
    template = read_json(ROOT / "app/resources/examples/rotating-shaft.kutta.json")
    template["mesh"] = mesh.name
    for speed in SPEEDS:
        template["boundaries"]["shaft"]["angular_speed"] = speed
        write_json(root / f"shaft-{speed}.kutta.json", template)

    subprocess.run([cli, "solve", root / "shaft-8.kutta.json", "--backend", "cpu",
                    "--iterations", "1", "--out", root / "prime"],
                   stdout=subprocess.DEVNULL, check=True)
    groups = {"cold": [], "forward": [], "reverse": []}
    # Interleave the chains so a drifting host load is shared across the modes.
    for i, speed in enumerate(SPEEDS):
        for mode, omega in (("cold", speed), ("forward", speed), ("reverse", SPEEDS[-1-i])):
            runs = groups[mode]
            out = root / f"{mode}-{omega}"
            args = [cli, "solve", root / f"shaft-{omega}.kutta.json", "--backend", "cpu",
                    "--progress", "jsonl", "--iterations", "4000", "--write-restart", "--out", out]
            restart = mode != "cold" and bool(runs)
            if restart:
                args.extend(["--restart", root / f"{mode}-{runs[-1]['omega']}" / "restart.bin"])
            with (root / f"{out.name}.jsonl").open("w") as log:
                before = time.perf_counter()
                result = subprocess.run(args, stdout=log, stderr=subprocess.PIPE, text=True,
                                        env={**os.environ, "KUTTA_MEMORY_PHASES": "1"}, timeout=120)
                elapsed = time.perf_counter() - before
            (root / f"{out.name}.stderr").write_text(result.stderr)
            result.check_returncode()
            summary = read_json(out / "summary.json")
            if summary["verdict"] != "converged" or summary["convergence"]["residuals_waived"]:
                raise RuntimeError(f"{out.name}: {summary['verdict']}")
            initialization = summary["initialization"]
            if initialization != ("native-restart" if restart else "standard"):
                raise RuntimeError(f"{out.name}: unexpected initialization {initialization}")
            peaks = [int(n) for n in re.findall(r"^kutta-memory .*? peak=(\d+)", result.stderr, re.M)]
            if not peaks or not max(peaks):
                raise RuntimeError(f"{out.name}: missing process peak memory")
            loads = read_json(out / "loads.json")["zones"]
            runs.append({"omega": omega, "iterations": summary["iterations"], "elapsed_s": elapsed,
                         "solver_s": summary["wall_time_seconds"], "peak_process_bytes": max(peaks),
                         "shaft_torque_nm": torque(loads, "shaft"), "housing_torque_nm": torque(loads, "housing"),
                         "restart_source_sha256": summary["restart_source_sha256"],
                         "verdict": summary["verdict"]})
            print(out.name, runs[-1], flush=True)

    cold = {run["omega"]: run for run in groups["cold"]}
    for mode in ("forward", "reverse"):
        for run in groups[mode]:
            baseline = cold[run["omega"]]
            for field in ("shaft_torque_nm", "housing_torque_nm"):
                run[field + "_delta_vs_cold"] = run[field] - baseline[field]
    write_json(root / "sweep.json", {"host": platform.platform(), "backend": "cpu",
                                    "solver_sha256": sha256(cli), "mesh_sha256": sha256(mesh),
                                    "termination": template["monitors"], "iteration_cap": 4000,
                                    "runs": groups})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("solver", type=Path)
    parser.add_argument("output", type=Path, help="fresh directory for all fifteen runs")
    args = parser.parse_args()
    sweep(args.solver.resolve(), args.output.resolve())
