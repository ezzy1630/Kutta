"""Compare decaying and sustained SST freestream turbulence at two Zerihan ride heights."""

import argparse
import csv
import hashlib
import json
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def run(solver, out):
    template = ROOT / "validation/cases/s1223-uiuc.kutta.json"
    data = (
        ROOT
        / "validation/data/zerihan_tyrrell026_ground_effect/table15_2d_sectional_Cl_vs_hc.csv"
    )
    with data.open() as src:
        observed = {
            r["h_over_c"]: float(r["Cl_expt_centre_pressure_integrated"])
            for r in csv.DictReader(src)
        }
    results = []
    for height in ("0.224", "0.090"):
        for sustain in (False, True):
            case = json.loads(template.read_text())
            case["mesh"] = str(ROOT / f"runs/fsae_aero/mesh/zerihan_h{height}_L2.msh")
            case["material"]["viscosity"] = case["reference_values"]["viscosity"] = (
                1 / 462000
            )
            for name in ("inlet", "outlet"):
                case["boundaries"][name]["turbulence"]["intensity"] = 0.002
            case["boundaries"]["ground"] = {
                "type": "wall",
                "motion": "moving",
                "direction": [1, 0, 0],
                "speed": 1,
            }
            case["reports"]["lift"]["direction"] = [0, -1, 0]
            k = 1.5 * 0.002**2
            ambient = {"k": k, "omega": k * 462000}
            case["initialization"].update(ambient)
            if sustain:
                case["models"]["k_omega_sst"] = {"ambient": ambient}
            d = out / f"h{height}-{'sustained' if sustain else 'decaying'}"
            d.mkdir(parents=True, exist_ok=True)
            path = d / "case.kutta.json"
            path.write_text(json.dumps(case, indent=2))
            with (d / "run.log").open("w") as log:
                code = subprocess.run(
                    [
                        solver,
                        "solve",
                        path,
                        "--backend",
                        "metal",
                        "--iterations",
                        "6000",
                        "--out",
                        d,
                    ],
                    check=False,
                    stdout=log,
                    stderr=subprocess.STDOUT,
                ).returncode
            item = {
                "h_over_c": float(height),
                "sustained": sustain,
                "measured_sectional_estimate": observed[height],
                "exit_code": code,
                "solver_sha256": hashlib.sha256(solver.read_bytes()).hexdigest(),
            }
            if not code:
                summary = json.loads((d / "summary.json").read_text())
                item.update(
                    {
                        k: summary[k]
                        for k in ("mesh_sha256", "case_sha256", "iterations", "verdict")
                    }
                )
                with (d / "forces.csv").open() as src:
                    item["coefficients"] = {
                        r["name"]: float(r["coefficient"])
                        for r in csv.DictReader(src)
                        if r["kind"] == "report"
                    }
            results.append(item)
            (out / "comparison.json").write_text(json.dumps(results, indent=2) + "\n")
            print(item, flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("solver", type=Path)
    p.add_argument("out", type=Path)
    a = p.parse_args()
    run(a.solver.resolve(), a.out.resolve())
