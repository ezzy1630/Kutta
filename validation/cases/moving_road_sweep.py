"""Compare cold and native warm starts on the same two-element wing and moving road."""

import argparse
import csv
import os
import platform
import re
import subprocess
import time
from pathlib import Path

from native_restart_sweep import read_json, sha256, write_json


def compare(cli, setup, root, backend):
    root.mkdir(parents=True)
    case = read_json(setup)
    mesh = (setup.parent / case["mesh"]).resolve()
    case["mesh"] = str(mesh)
    binary_hash, mesh_hash = sha256(cli), sha256(mesh)
    inlet, road = (case["boundaries"][name] for name in ("inlet", "ground"))
    speed = inlet["speed"]
    if (
        road.get("motion") != "moving"
        or road.get("speed") != speed
        or case["reference_values"]["velocity"] != speed
    ):
        raise ValueError("The inlet, moving ground and reference speed must match")
    speeds = [round(speed * ratio, 8) for ratio in (0.9, 0.95, 1.0, 1.05, 1.1)]
    for value in speeds:
        inlet["speed"] = road["speed"] = case["reference_values"]["velocity"] = value
        write_json(root / f"speed-{value}.kutta.json", case)
    groups = {mode: [] for mode in ("cold", "forward", "reverse")}
    for i, speed in enumerate(speeds):
        for mode, value in (
            ("cold", speed),
            ("forward", speed),
            ("reverse", speeds[-1 - i]),
        ):
            out = root / f"{mode}-{value}"
            cmd = [
                cli,
                "solve",
                root / f"speed-{value}.kutta.json",
                "--backend",
                backend,
                "--iterations",
                "1600",
                "--write-restart",
                "--out",
                out,
            ]
            warm = mode != "cold" and bool(groups[mode])
            if warm:
                cmd.extend(
                    [
                        "--restart",
                        root
                        / f"{mode}-{groups[mode][-1]['speed_m_s']}"
                        / "restart.bin",
                    ]
                )
            if sha256(cli) != binary_hash:
                raise RuntimeError("Solver executable changed during the sweep")
            started = time.perf_counter()
            with (root / f"{mode}-{value}.log").open("w") as log:
                result = subprocess.run(
                    cmd,
                    check=False,
                    stdout=log,
                    stderr=subprocess.PIPE,
                    text=True,
                    env={**os.environ, "KUTTA_MEMORY_PHASES": "1"},
                )
            elapsed = time.perf_counter() - started
            (root / f"{mode}-{value}.stderr").write_text(result.stderr)
            result.check_returncode()
            summary = read_json(out / "summary.json")
            if summary["mesh_sha256"] != mesh_hash or warm != (
                summary["initialization"] == "native-restart"
            ):
                raise ValueError("Unexpected mesh or initialization")
            with (out / "forces.csv").open() as source:
                loads = {
                    r["name"]: float(r["total_projection"])
                    for r in csv.DictReader(source)
                    if r["kind"] == "report"
                }
            if set(loads) != {"drag", "downforce"}:
                raise ValueError("Expected dimensional drag and downforce reports")
            peaks = re.findall(r"^kutta-memory .*? peak=(\d+)", result.stderr, re.MULTILINE)
            groups[mode].append(
                {
                    "speed_m_s": value,
                    "elapsed_s": elapsed,
                    "loads_N": loads,
                    "peak_process_bytes": max(map(int, peaks)),
                    **{
                        k: summary[k]
                        for k in (
                            "initialization",
                            "iterations",
                            "verdict",
                            "case_sha256",
                            "restart_source_sha256",
                        )
                    },
                    "residuals_waived": summary["convergence"]["residuals_waived"],
                }
            )
            write_json(
                root / "sweep.json",
                {
                    "mesh_sha256": mesh_hash,
                    "solver_sha256": binary_hash,
                    "backend": backend,
                    "platform": platform.platform(),
                    "groups": groups,
                },
            )
            print(mode, value, summary["verdict"], summary["iterations"], flush=True)
    if any(
        x["verdict"] != "converged" or x["residuals_waived"]
        for group in groups.values()
        for x in group
    ):
        raise RuntimeError("Some operating points missed the common stopping criteria")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("cli", "setup", "out"):
        parser.add_argument(name, type=Path)
    parser.add_argument("--backend", default="cpu")
    args = parser.parse_args()
    compare(args.cli.resolve(), args.setup.resolve(), args.out.resolve(), args.backend)
