"""Physical time refinement for the impulsively started concentric shaft, using the cylindrical diffusion eigenmodes."""

import argparse
import csv
import json
import math
import subprocess
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import j1, y1

ROOT = Path(__file__).resolve().parents[2]


def exact(radius, time):
    ri, ro, nu, omega = 0.1, 0.2, 0.01, 10.0

    def steady(r):
        return omega * ri**2 / (ro**2 - ri**2) * (ro**2 / r - r)

    def mode(lam, r):
        return j1(lam * r) * y1(lam * ri) - y1(lam * r) * j1(lam * ri)

    result = steady(radius)
    for m in range(1, 31):
        center = m * math.pi / (ro - ri)
        half = 0.45 * math.pi / (ro - ri)
        lam = brentq(lambda k: mode(k, ro), center - half, center + half)
        numerator = quad(lambda r, lam=lam: r * steady(r) * mode(lam, r), ri, ro)[0]
        denominator = quad(lambda r, lam=lam: r * mode(lam, r) ** 2, ri, ro)[0]
        result -= (
            (numerator / denominator)
            * mode(lam, radius)
            * math.exp(-nu * lam**2 * time)
        )
    return result


def run(cli, out, backend):
    case = ROOT / "app/resources/examples/rotating-shaft.kutta.json"
    out.mkdir(parents=True, exist_ok=True)
    post = out / "shaft-line.json"
    post.write_text(
        json.dumps(
            {
                "v": 1,
                "probes": [
                    {
                        "type": "line",
                        "start": [0.102, 0, 0],
                        "end": [0.198, 0, 0],
                        "count": 65,
                        "fields": ["velocity_y"],
                    }
                ],
            }
        )
    )
    values, results = [], []
    for steps in (4, 8, 16, 32):
        dt = 0.04 / steps
        directory = out / f"n{steps}"
        subprocess.run(
            [
                cli,
                "transient",
                case,
                "--backend",
                backend,
                "--dt",
                str(dt),
                "--steps",
                str(steps),
                "--tolerance",
                "0.000003",
                "--correctors",
                "160",
                "--out",
                directory,
            ],
            check=True,
        )
        subprocess.run([cli, "export-view", directory, "--post", post], check=True)
        with (directory / "view/probe-0.csv").open() as src:
            data = list(csv.DictReader(src))
        if any(r["valid"] != "1" for r in data):
            raise ValueError("Shaft samples outside the annulus")
        u = np.array([float(r["velocity_y"]) for r in data])
        ref = exact(np.array([float(r["x"]) for r in data]), steps * dt)
        rel = np.linalg.norm(u - ref) / np.linalg.norm(ref)
        values.append(u)
        results.append({"steps": steps, "dt_s": dt, "radial_profile_relative_l2": rel})
    for i in range(2, len(values)):
        ratio = np.linalg.norm(values[i - 2] - values[i - 1]) / np.linalg.norm(
            values[i - 1] - values[i]
        )
        results[i]["temporal_self_convergence_order"] = math.log2(ratio)
    (out / "comparison.json").write_text(json.dumps(results, indent=2) + "\n")
    print(results)
    if (
        results[-1]["radial_profile_relative_l2"] > 0.02
        or min(r["temporal_self_convergence_order"] for r in results[2:]) < 1.6
    ):
        raise SystemExit("Temporal or analytical profile accuracy gate missed")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("cli", type=Path)
    p.add_argument("out", type=Path)
    p.add_argument("--backend", default="cpu")
    a = p.parse_args()
    run(a.cli.resolve(), a.out.resolve(), a.backend)
