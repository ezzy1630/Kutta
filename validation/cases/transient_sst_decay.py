"""Verify the native transient SST fields, their sampled mean, and freestream sustaining terms."""

import argparse
import csv
import hashlib
import json
import math
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
POINTS = [[0.1, 0.1, 0.1], [0.9, 0.7, 0.5]]


def exact(t, sustained):
    factor = 1 if sustained else 1 + 0.0828 * 10 * t
    return {"k": 0.05 * factor ** (-0.09 / 0.0828), "omega": 10 / factor}


def profile_error(cli, run, spec, expected):
    subprocess.run([cli, "export-view", run, "--post", spec], check=True)
    with (run / "view/probe-0.csv").open() as source:
        rows = list(csv.DictReader(source))
    if len(rows) != len(POINTS) or any(r["valid"] != "1" for r in rows):
        raise ValueError("Missing SST verification probes")
    return {
        key: max(abs(float(r[key]) / value - 1) for r in rows)
        for key, value in expected.items()
    }


def run(cli, root, backend):
    root.mkdir(parents=True, exist_ok=True)
    mesh = ROOT / "tests/data/meshes/fsae_box.msh"
    spec = root / "probe.json"
    spec.write_text(
        json.dumps(
            {
                "v": 1,
                "probes": [
                    {"type": "points", "points": POINTS, "fields": ["k", "omega"]}
                ],
            }
        )
    )
    boundaries = {
        name: {"type": "symmetry"}
        for name in ("inlet", "front_wing", "ground", "top", "side_left", "side_right")
    }
    boundaries["outlet"] = {"type": "pressure-outlet", "gauge_pressure": 0}
    results = []
    for sustained, steps in [(False, 4), (False, 8), (False, 16), (True, 16)]:
        case = {
            "kutta_case": 1,
            "mesh": str(mesh),
            "models": {"viscous": "k-omega-sst"},
            "material": {"density": 1, "viscosity": 1e-5},
            "boundaries": boundaries,
            "initialization": {
                "method": "standard",
                "velocity": [0, 0, 0],
                **exact(0, False),
            },
            "methods": {"pressure_velocity_coupling": "simplec"},
        }
        if sustained:
            case["models"]["k_omega_sst"] = {"ambient": exact(0, True)}
        tag = f"{'sustained' if sustained else 'decay'}-{steps}"
        setup, out = root / (tag + ".kutta.json"), root / tag
        setup.write_text(json.dumps(case))
        dt = 0.4 / steps
        stride = 2 if steps == 8 else 1
        subprocess.run(
            [
                cli,
                "transient",
                setup,
                "--backend",
                backend,
                "--dt",
                str(dt),
                "--steps",
                str(steps),
                "--correctors",
                "120",
                "--tolerance",
                "2e-7",
                "--mean-every",
                str(stride),
                "--mean-from",
                "1",
                "--out",
                out,
            ],
            check=True,
        )
        final = profile_error(cli, out, spec, exact(0.4, sustained))
        sampled = [exact(j * dt, sustained) for j in range(1, steps + 1, stride)]
        summary = json.loads((out / "mean/summary.json").read_text())
        if summary["samples"] != len(sampled) or not math.isclose(
            summary["last_sample_time_s"], (1 + (len(sampled) - 1) * stride) * dt
        ):
            raise ValueError("Mean-field sampling window differs from the request")
        mean = profile_error(
            cli,
            out / "mean",
            spec,
            {k: sum(row[k] for row in sampled) / len(sampled) for k in final},
        )
        results.append(
            {
                "sustained": sustained,
                "dt_s": dt,
                "sample_stride_steps": stride,
                "final_relative_error": final,
                "mean_relative_error": mean,
            }
        )
    orders = {
        k: [
            math.log2(
                results[i - 1]["final_relative_error"][k]
                / results[i]["final_relative_error"][k]
            )
            for i in (1, 2)
        ]
        for k in ("k", "omega")
    }
    record = {
        "solver_sha256": hashlib.sha256(cli.read_bytes()).hexdigest(),
        "backend": backend,
        "mesh_sha256": hashlib.sha256(mesh.read_bytes()).hexdigest(),
        "samples": results,
        "decay_orders": orders,
    }
    (root / "comparison.json").write_text(json.dumps(record, indent=2) + "\n")
    if (
        min(min(v) for v in orders.values()) < 1.6
        or max(results[2]["mean_relative_error"].values()) > 0.001
        or max(results[3]["final_relative_error"].values()) > 1e-5
    ):
        raise ValueError(
            "SST decay, mean-field accuracy or freestream equilibrium check failed"
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cli", type=Path)
    parser.add_argument("out", type=Path, nargs="?")
    parser.add_argument("--backend", default="cpu")
    args = parser.parse_args()
    if args.out:
        run(args.cli.resolve(), args.out.resolve(), args.backend)
    else:
        with tempfile.TemporaryDirectory() as directory:
            run(args.cli.resolve(), Path(directory), args.backend)
