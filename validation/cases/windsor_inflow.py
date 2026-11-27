"""Sample the Windsor inlet approach flow and integrate the local floor boundary-layer profile."""

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


def sample(cli, run):
    summary = json.loads((run / "summary.json").read_text())
    if (
        summary["mesh_sha256"]
        != "c8c33dc8c6d07af42de1843e56095ec7710613cce6fc345ced5421839ebb6058"
    ):
        raise ValueError("Expected the public AutoCFD5 G1 mesh")
    with tempfile.TemporaryDirectory() as d:
        root = Path(d)
        for name in ("summary.json", "view-input.bin"):
            try:
                os.link(run / name, root / name)
            except OSError:
                shutil.copyfile(run / name, root / name)
        fields = ["velocity_x", "velocity_magnitude", "k", "omega"]
        probes = [
            {
                "name": "Approach profile",
                "type": "line",
                "start": [-1, 0, 0],
                "end": [-1, 0, 0.5],
                "count": 501,
                "fields": fields,
            },
            {
                "name": "Freestream",
                "type": "points",
                "points": [[-2, 0, 1.3], [-1, 0, 0.5], [0, 0, 1.15]],
                "fields": fields,
            },
        ]
        spec = root / "post.json"
        spec.write_text(json.dumps({"v": 1, "probes": probes}))
        subprocess.run([cli, "export-view", root, "--post", spec], check=True)

        def read(i):
            with (root / "view" / f"probe-{i}.csv").open() as src:
                return list(csv.DictReader(src))

        line, free = read(0), read(1)
    if any(r["valid"] != "1" for r in line + free):
        raise ValueError("An approach-flow probe is outside the volume")
    edge = float(free[1]["velocity_x"])
    points = [(float(r["z"]), float(r["velocity_x"]) / edge) for r in line]
    # Impose the known no-slip value at the floor if cell reconstruction differs.
    points[0] = (0.0, 0.0)
    for i, ((z0, r0), (z1, r1)) in enumerate(itertools.pairwise(points)):
        if r0 < 0.99 <= r1:
            z99 = z0 + (0.99 - r0) * (z1 - z0) / (r1 - r0)
            profile = points[: i + 1] + [(z99, 0.99)]
            break
    else:
        raise ValueError("The 99% edge is outside the sampled profile")

    def integrate(f):
        return sum(
            (z1 - z0) * (f(r0) + f(r1)) / 2
            for (z0, r0), (z1, r1) in itertools.pairwise(profile)
        )

    return {
        "mesh_sha256": summary["mesh_sha256"],
        "case_sha256": summary["case_sha256"],
        "approach_x_m": -1,
        "edge_reference_z_m": 0.5,
        "edge_reference_velocity_m_s": edge,
        "velocity_ratio_by_z_m": {
            str(points[i][0]): points[i][1] for i in (10, 30, 60, 100, 200, 300, 500)
        },
        "delta_99_m": z99,
        "displacement_to_delta99_m": integrate(lambda r: 1 - r),
        "momentum_to_delta99_m": integrate(lambda r: r * (1 - r)),
        "isotropic_Tu_probes": [
            {
                "xyz_m": [float(r[c]) for c in "xyz"],
                "velocity_m_s": float(r["velocity_magnitude"]),
                "k": float(r["k"]),
                "omega": float(r["omega"]),
                "Tu": math.sqrt(2 * float(r["k"]) / 3) / float(r["velocity_magnitude"]),
            }
            for r in free
        ],
    }


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("cli", type=Path)
    p.add_argument("run", type=Path)
    p.add_argument("out", type=Path)
    a = p.parse_args()
    a.out.write_text(
        json.dumps(sample(a.cli.resolve(), a.run.resolve()), indent=2) + "\n"
    )
