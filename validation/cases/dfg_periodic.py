"""Generate and run the FeatFlow DFG 2D-2 periodic cylinder benchmark (SI units)."""

import argparse
import csv
import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np
from scipy.signal import find_peaks

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools/fsae_aero"))
import itertools

import fluent_msh


def generate(root, scale):
    h, r = 0.41, 0.05
    side = max(8, round(28 / scale))
    radial = max(6, round(16 / scale))
    nx = max(12, round(1.79 / (0.018 * scale)))
    corners = np.array([[0, 0], [h, 0], [h, h], [0, h]])
    points = []
    cells = []
    channel = []
    cylinder = []
    inlet = []
    outlet = []
    rings = []
    for j in range(radial + 1):
        ring = []
        t = (math.exp(2 * j / radial) - 1) / math.expm1(2)
        for sector in range(4):
            for k in range(side):
                q = k / side
                theta = math.pi * (-0.75 + (sector + q) * 0.5)
                inner = np.array([0.2 + r * math.cos(theta), 0.2 + r * math.sin(theta)])
                outer = (1 - q) * corners[sector] + q * corners[(sector + 1) % 4]
                ring.append(len(points))
                points.append([*((1 - t) * inner + t * outer), 0])
        rings.append(ring)
    for j in range(radial):
        for i in range(4 * side):
            k = (i + 1) % (4 * side)
            cells.append([rings[j][i], rings[j + 1][i], rings[j + 1][k], rings[j][k]])
    for i in range(4 * side):
        k = (i + 1) % (4 * side)
        cylinder.append([rings[0][i], rings[0][k]])
        edge = [rings[-1][i], rings[-1][k]]
        if i // side in (0, 2):
            channel.append(edge)
        if i // side == 3:
            inlet.append(edge)
    column = rings[-1][side : 2 * side + 1]
    for ix in range(1, nx + 1):
        new = []
        for k in range(side + 1):
            new.append(len(points))
            points.append([h + (2.2 - h) * ix / nx, h * k / side, 0])
        for k in range(side):
            cells.append([column[k], new[k], new[k + 1], column[k + 1]])
        channel.extend([[column[0], new[0]], [column[-1], new[-1]]])
        column = new
    outlet = list(itertools.pairwise(column))
    zones = [
        ("channel", "wall", np.array(channel)),
        ("cylinder", "wall", np.array(cylinder)),
        ("outlet", "pressure-outlet", np.array(outlet)),
    ]
    bc = {
        "channel": {"type": "wall"},
        "cylinder": {"type": "wall"},
        "outlet": {
            "type": "pressure-outlet",
            "gauge_pressure": 0,
            "backflow_pressure": "static",
        },
    }
    for i, e in enumerate(inlet):
        name = f"inlet-{i}"
        zones.append((name, "velocity-inlet", np.array([e])))
        y0, y1 = points[e[0]][1], points[e[1]][1]
        speed = 6 * ((y0 + y1) * h / 2 - (y0 * y0 + y0 * y1 + y1 * y1) / 3) / h**2
        bc[name] = {"type": "velocity-inlet", "velocity": [speed, 0, 0]}
    fluent_msh.write(
        root / "cylinder.msh", np.array(points), [(3, np.array(cells))], zones
    )
    case = {
        "kutta_case": 1,
        "mesh": "cylinder.msh",
        "models": {"viscous": "laminar"},
        "material": {"density": 1, "viscosity": 0.001},
        "boundaries": bc,
        "reference_values": {
            "length": 0.1,
            "area": 0.1,
            "density": 1,
            "velocity": 1,
            "viscosity": 0.001,
        },
        "reports": {
            name: {
                "type": name,
                "zones": ["cylinder"],
                "direction": axis,
                "coefficient": True,
            }
            for name, axis in (("drag", [1, 0, 0]), ("lift", [0, 1, 0]))
        },
        "methods": {
            "pressure_velocity_coupling": "simplec",
            "gradient_limiter": "differentiable",
        },
        "initialization": {"method": "standard", "velocity": [0, 0, 0]},
        "controls": {"relaxation": {"momentum": 0.9, "pressure": 1}},
    }
    path = root / "case.kutta.json"
    path.write_text(json.dumps(case, indent=2))
    return path


def summarize(path):
    with (path / "physical-time.csv").open() as src:
        rows = list(csv.DictReader(src))
    t = np.array([float(r["time_s"]) for r in rows])
    cl = np.array([float(r["lift [N]"]) / 0.05 for r in rows])
    cd = np.array([float(r["drag [N]"]) / 0.05 for r in rows])
    peaks = find_peaks(cl, prominence=0.15)[0]
    if len(peaks) < 5:
        raise ValueError("Fewer than four complete shedding cycles")
    selected = peaks[-5:]
    sl = slice(selected[0], selected[-1] + 1)
    period = np.diff(t[selected]).mean()
    result = {
        "cycles": 4,
        "from_s": t[selected[0]],
        "to_s": t[selected[-1]],
        "St": 0.1 / period,
        "max_CL": cl[sl].max(),
        "max_CD": cd[sl].max(),
        "peak_lift_by_cycle": cl[selected].tolist(),
    }
    if t[-1] >= 30 - 1e-9:
        troughs = find_peaks(-cl, prominence=0.15)[0]
        troughs = troughs[(t[troughs] >= 25) & (t[troughs] <= 30)]
        if len(troughs) < 2:
            raise ValueError(
                "No complete lift-minimum cycle in the 25–30 s reference window"
            )
        start, end = troughs[-2:]
        result["reference_window_cycle"] = {
            "from_s": t[start],
            "to_s": t[end],
            "St": 0.1 / (t[end] - t[start]),
        }
        for name, series in (("CD", cd), ("CL", cl)):
            lo, hi = series[start : end + 1].min(), series[start : end + 1].max()
            result["reference_window_cycle"][name] = {
                "min": lo,
                "max": hi,
                "midrange": (lo + hi) / 2,
                "amplitude": hi - lo,
            }
    return result


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("cli", type=Path)
    p.add_argument("out", type=Path)
    p.add_argument("--scale", type=float, default=1)
    p.add_argument("--dt", type=float, default=0.005)
    p.add_argument("--steps", type=int, help="Override the default 30 s integration")
    p.add_argument("--backend", default="cpu")
    p.add_argument("--correctors", type=int, default=90)
    p.add_argument(
        "--steady-seed-iterations",
        type=int,
        default=0,
        help="Initialize from this many steady iterations on the same mesh",
    )
    a = p.parse_args()
    if a.dt <= 0 or a.scale <= 0:
        p.error("--dt and --scale must be positive")
    if a.steady_seed_iterations < 0:
        p.error("--steady-seed-iterations must be nonnegative")
    if a.steps is None:
        a.steps = math.ceil(30 / a.dt)
    a.out.mkdir(parents=True, exist_ok=True)
    case = generate(a.out, a.scale)
    out = a.out / "result"
    restart = []
    if a.steady_seed_iterations:
        start = a.out / "start"
        subprocess.run(
            [
                a.cli,
                "solve",
                case,
                "--backend",
                a.backend,
                "--iterations",
                str(a.steady_seed_iterations),
                "--write-restart",
                "--out",
                start,
            ],
            check=True,
        )
        restart = ["--restart", start / "restart.bin"]
    subprocess.run(
        [
            a.cli,
            "transient",
            case,
            "--dt",
            str(a.dt),
            "--steps",
            str(a.steps),
            "--correctors",
            str(a.correctors),
            "--tolerance",
            "0.0001",
            "--backend",
            a.backend,
            *restart,
            "--out",
            out,
        ],
        check=True,
    )
    (out / "shedding.json").write_text(json.dumps(summarize(out), indent=2) + "\n")
