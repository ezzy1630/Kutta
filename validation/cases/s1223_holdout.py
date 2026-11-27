"""Run unchanged Kutta wing methods at two Reynolds numbers outside the earlier polar study."""

import argparse
import csv
import hashlib
import itertools
import json
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ANGLE = 4.76


def reference(reynolds):
    with (
        ROOT / "validation/data/uiuc_lsat_s1223/s1223_vol1_drag_polars.csv"
    ).open() as src:
        rows = sorted(
            (float(r["alpha_deg"]), float(r["cl"]), float(r["cd"]))
            for r in csv.DictReader(src)
            if r["config"] == "clean" and int(r["re_avg"]) == reynolds
        )
    for low, high in itertools.pairwise(rows):
        if low[0] <= ANGLE <= high[0]:
            t = (ANGLE - low[0]) / (high[0] - low[0])
            return {
                k: low[i] + t * (high[i] - low[i])
                for i, k in enumerate(("alpha", "lift", "drag"))
            }
    raise ValueError("Angle is not bracketed by this measured polar")


def run(solver, out, backend, models):
    template = ROOT / "validation/cases/s1223-uiuc.kutta.json"
    binary = hashlib.sha256(solver.read_bytes()).hexdigest()
    results = []
    for reynolds in (149700, 249700):
        for model in models:
            case = json.loads(template.read_text())
            case["mesh"] = str((template.parent / case["mesh"]).resolve())
            case["material"]["viscosity"] = case["reference_values"]["viscosity"] = (
                1 / reynolds
            )
            case["initialization"]["omega"] = case["initialization"]["k"] * reynolds
            if model == "sustained":
                case["models"]["k_omega_sst"] = {
                    "ambient": {
                        key: case["initialization"][key] for key in ("k", "omega")
                    }
                }
            if model == "intermittency":
                case["models"]["transition"] = model
            directory = out / f"re{reynolds}-{model}"
            directory.mkdir(parents=True, exist_ok=True)
            path = directory / "case.kutta.json"
            path.write_text(json.dumps(case, indent=2))
            with (directory / "solver.log").open("w") as log:
                status = subprocess.run(
                    [
                        solver,
                        "solve",
                        path,
                        "--backend",
                        backend,
                        "--iterations",
                        "6000",
                        "--out",
                        directory,
                    ],
                    check=False,
                    stdout=log,
                    stderr=subprocess.STDOUT,
                )
            result = {
                "Re": reynolds,
                "alpha_deg": ANGLE,
                "model": model,
                "binary_sha256": binary,
                "measured_interpolated": reference(reynolds),
                "exit_status": status.returncode,
            }
            if status.returncode == 0:
                summary = json.loads((directory / "summary.json").read_text())
                with (directory / "forces.csv").open() as src:
                    result["predicted"] = {
                        r["name"]: float(r["coefficient"])
                        for r in csv.DictReader(src)
                        if r["kind"] == "report"
                    }
                for key in (
                    "verdict",
                    "iterations",
                    "mesh_sha256",
                    "case_sha256",
                    "wall_time_seconds",
                ):
                    result[key] = summary[key]
            results.append(result)
            (out / "comparison.json").write_text(json.dumps(results, indent=2) + "\n")
            print(reynolds, model, result.get("verdict", "failed"), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("solver", type=Path)
    parser.add_argument("out", type=Path)
    parser.add_argument("--backend", default="metal")
    parser.add_argument(
        "--models",
        nargs="+",
        choices=["sst", "intermittency", "sustained"],
        default=["sst", "intermittency", "sustained"],
    )
    args = parser.parse_args()
    run(args.solver.resolve(), args.out.resolve(), args.backend, args.models)
