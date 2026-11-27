"""Cross-check a 2-D S1223 wall drag against a simulated far-wake total-head traverse."""

import argparse
import csv
import itertools
import json
import math
import os
import shutil
import subprocess
import tempfile
from pathlib import Path


def check(cli, run):
    if (
        json.loads((run / "summary.json").read_text())["mesh_sha256"]
        != "d762ed53700668caf659dc8cd785705af42182b0b838e77456b06328c613e034"
    ):
        raise ValueError("Expected the unit-chord S1223 mesh at 4.76 degrees")
    ref = 1.0  # unit chord, U and rho in the native S1223 case
    q = 0.5
    alpha = math.radians(4.76)
    # The native geometry rotates about quarter chord.
    x_te = 0.25 + 0.75 * math.cos(alpha)
    requests = []
    for dx in (1.0, 1.25, 1.5):
        x = x_te + dx
        requests.append(
            {
                "type": "line",
                "name": f"Wake at TE + {dx} chord",
                "start": [x, -0.65, 0],
                "end": [x, 0.35, 0],
                "count": 601,
                "fields": ["total_pressure", "pressure", "velocity_x", "velocity_y"],
            }
        )
    requests.append(
        {
            "type": "line",
            "name": "Inlet reference",
            "start": [-29, 0, 0],
            "end": [-29, 0.01, 0],
            "count": 2,
            "fields": ["total_pressure", "pressure", "velocity_x"],
        }
    )
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        for name in ("summary.json", "view-input.bin"):
            try:
                os.link(run / name, root / name)
            except OSError:
                shutil.copyfile(run / name, root / name)
        post = root / "post.json"
        post.write_text(json.dumps({"v": 1, "probes": requests}))
        subprocess.run([cli, "export-view", root, "--post", post], check=True)

        def read(i):
            with (root / "view" / f"probe-{i}.csv").open() as src:
                data = list(csv.DictReader(src))
            if any(r["valid"] != "1" for r in data):
                raise ValueError("Wake traverse falls outside the mesh")
            return data

        inlet = read(3)[0]
        pt = float(inlet["total_pressure"])
        results = []
        for i, dx in enumerate((1.0, 1.25, 1.5)):
            rows = read(i)
            loss = [(pt - float(r["total_pressure"])) / q for r in rows]
            if max(loss) > 1:
                raise ValueError(
                    "Far-wake inversion invalid: total head below freestream static"
                )
            speed = [math.sqrt(1 - d) for d in loss]
            integrand = [2 * u * (1 - u) for u in speed]
            y = [float(r["y"]) for r in rows]
            cd = sum(
                (y[j + 1] - y[j]) * (a + b) * 0.5 / ref
                for j, (a, b) in enumerate(itertools.pairwise(integrand))
            )
            results.append(
                {
                    "x_minus_te_over_c": dx,
                    "wake_head_loss_CD": cd,
                    "outer_edge_loss_over_q": [loss[0], loss[-1]],
                    "max_static_offset_over_q": max(
                        abs(float(r["pressure"]) - float(inlet["pressure"])) / q
                        for r in rows
                    ),
                }
            )
    with (run / "forces.csv").open() as src:
        drag = next(
            r
            for r in csv.DictReader(src)
            if r["kind"] == "report" and r["name"] == "drag"
        )
    return {
        "wall_CD": float(drag["coefficient"]),
        "wall_pressure_CD": float(drag["pressure_projection"]) / q,
        "wall_viscous_CD": float(drag["viscous_projection"]) / q,
        "inlet_total_pressure_Pa": pt,
        "traverses": results,
    }


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("cli", type=Path)
    p.add_argument("runs", type=Path, nargs="+")
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    a.out.write_text(
        json.dumps({d.name: check(a.cli.resolve(), d) for d in a.runs}, indent=2) + "\n"
    )
