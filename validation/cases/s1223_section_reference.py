"""Independent viscous XFOIL cross-check on the UIUC as-tested S1223 section."""

import argparse
import csv
import hashlib
import json
import math
import subprocess
import tempfile
from pathlib import Path


def evaluate(executable, coordinates, reynolds, angle):
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        (root / "section.dat").write_bytes(coordinates.read_bytes())
        commands = f"PLOP\nG\n\nLOAD section.dat\nPANE\nOPER\nVISC {reynolds}\nITER 250\nPACC\npolar.txt\n\nALFA {angle}\nPACC\n\nQUIT\n"
        result = subprocess.run(
            [executable],
            cwd=root,
            input=commands,
            text=True,
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            timeout=60,
        )
        if result.returncode:
            raise RuntimeError(
                f"XFOIL failed at Re={reynolds}, incidence={angle}: {result.stdout[-400:]}"
            )
        rows = (root / "polar.txt").read_text().splitlines()
        separator = next(
            i for i, line in enumerate(rows) if line.strip().startswith("-----")
        )
        data = [line.split() for line in rows[separator + 1 :] if line.strip()]
        if not data:
            return None
        a, cl, cd, _, _, xtop, xbot = map(float, data[0][:7])
        if (
            len(data) != 1
            or abs(a - angle) > 0.001
            or not all(map(math.isfinite, (cl, cd, xtop, xbot)))
        ):
            raise ValueError("Unexpected XFOIL polar row")
        return {
            "lift": cl,
            "drag": cd,
            "transition_x_over_c_upper": xtop,
            "transition_x_over_c_lower": xbot,
        }


def run(executable, data, out):
    with (data / "s1223_vol1_drag_polars.csv").open() as source:
        rows = [
            r
            for r in csv.DictReader(source)
            if r["config"] == "clean"
            and min(
                abs(int(r["re_avg"]) - re) for re in (150000, 200000, 250000, 300000)
            )
            < 2500
        ]
    results = []
    coords = data / "lsat_vol1_COORD01_s1223.dap.txt"
    for r in rows:
        re, angle = int(r["re_avg"]), float(r["alpha_deg"])
        results.append(
            {
                "Re": re,
                "alpha_deg": angle,
                "measured_lift": float(r["cl"]),
                "measured_drag": float(r["cd"]),
                "predicted": evaluate(executable, coords, re, angle),
            }
        )
    evidence = {
        "executable_sha256": hashlib.sha256(executable.read_bytes()).hexdigest(),
        "coordinates_sha256": hashlib.sha256(coords.read_bytes()).hexdigest(),
        "method": "XFOIL 6.99, default panel distribution and Ncrit=9, free transition, independent cold starts",
        "results": results,
    }
    out.write_text(json.dumps(evidence, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("executable", "data", "out"):
        p.add_argument(name, type=Path)
    a = p.parse_args()
    run(a.executable.resolve(), a.data.resolve(), a.out.resolve())
