"""Calculate AutoCFD5 Windsor G1 loads in the body axes using the tunnel reference probe."""
import argparse
import csv
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from windsor_pressure import compare_sections, sections
from windsor_wake import compare_plane

# Loughborough Set B, Force_Mean/FCAV0D_NW_2_5yaw_Force.csv
MEASURED = {"drag": 0.3298, "lift": -0.0382, "side": 0.1509, "pitch": -0.0259}


def compare(cli, run, experiment=None, flow_data=None, wake_map=None):
    summary = json.loads((run / "summary.json").read_text())
    if summary["mesh_sha256"] != "c8c33dc8c6d07af42de1843e56095ec7710613cce6fc345ced5421839ebb6058":
        raise ValueError("Expected the AutoCFD5 Windsor G1 grid")
    probe = {"name": "Tunnel reference", "type": "line", "start": [-2, 0, 1.3], "end": [-1.99, 0, 1.3],
             "count": 2, "fields": ["pressure", "velocity_magnitude"]}
    with tempfile.TemporaryDirectory() as temp:
        workspace = Path(temp)
        for name in ("summary.json", "view-input.bin"):
            try:
                os.link(run / name, workspace / name)
            except OSError:
                shutil.copyfile(run / name, workspace / name)
        view = workspace / "view"
        spec = workspace / "post.json"

        def export(request):
            spec.write_text(json.dumps({"v": 1, **request}, separators=(",", ":")))
            subprocess.run([cli, "export-view", workspace, "--post", spec], check=True, stdout=sys.stderr)

        export({"probes": [probe]})
        with (view / "probe-0.csv").open() as source:
            reference = next(csv.DictReader(source))
        if reference["valid"] != "1":
            raise ValueError("Reference point falls outside the volume")
        speed = float(reference["velocity_magnitude"])
        pressure = float(reference["pressure"])
        yaw = math.radians(-2.5)
        drag = [math.cos(yaw), math.sin(yaw), 0]
        side = [-math.sin(yaw), math.cos(yaw), 0]
        walls = ["Flow.Windsor_Square_nW.Body", "Flow.Windsor_Square_nW.Base", "Flow.Windsor_Pins.Surface"]
        integral = {"name": "Windsor including pins", "zones": walls, "origin": [0, 0, 0],
                    "drag": drag, "side": side, "lift": [0, 0, 1], "fields": []}
        planes = sections(experiment, drag, side, walls) if experiment else []
        partitions = []
        for label, bounds in (("upper facing", (-1, -.5)), ("steep sides", (-.5, .5)), ("lower facing", (.5, 1))):
            partitions.append({**integral, "name": label, "zones": walls[:1],
                               "normal_range": {"axis": [0, 0, 1], "min": bounds[0], "max": bounds[1]}})
        export({"probes": [probe] + [section for section, _ in planes], "reference": {"velocity": speed, "pressure": pressure,
                "area": 0.112, "length": 0.6375}, "integrals": [integral] + partitions})
        taps = compare_sections(view, planes) if experiment else None
        report, *groups = json.loads((view / "reports.json").read_text())["reports"]
        wake = compare_plane(flow_data, drag, side, speed, export, view, wake_map) if flow_data else None
    total = report["total"]
    coefficients = total["coefficients"]["total"].copy()
    moment = total["moment_coefficients"]["total"]
    coefficients["pitch"] = sum(m * axis for m, axis in zip(moment, side))
    base = next(zone for zone in report["zones"] if zone["zone"] == "Flow.Windsor_Square_nW.Base")
    result = {"verdict": summary["verdict"], "iterations": summary["iterations"], "probe_speed_m_s": speed,
            "probe_pressure_Pa": pressure, "final_iterate_coefficients": coefficients,
            "base_pressure_drag_coefficient": base["coefficients"]["pressure"]["drag"],
            "measured_coefficients": MEASURED,
            "prediction_minus_measurement": {key: coefficients[key] - value for key, value in MEASURED.items()}}
    components = {group["name"]: group["total"]["coefficients"] for group in groups}
    components.update({zone["zone"]: zone["coefficients"] for zone in report["zones"] if zone["zone"] != walls[0]})
    for key in ("drag", "side", "lift"):
        if abs(sum(c["total"][key] for c in components.values()) - coefficients[key]) > 1e-8:
            raise ValueError(f"Surface partitions do not recover the total {key}")
    result["surface_load_coefficients"] = components
    if wake is not None:
        result["planar_wake"] = wake
    if taps is not None:
        result["surface_pressures"] = taps
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("solver", type=Path)
    parser.add_argument("run", type=Path)
    parser.add_argument("--pressure-data", type=Path, help="Extracted Pressure_Mean directory from Loughborough Set B")
    parser.add_argument("--flow-data", type=Path, help="Extracted FlowField_Mean directory from Loughborough Set B")
    parser.add_argument("--wake-map", type=Path, help="Save an experimental, CFD and difference plot (requires matplotlib)")
    args = parser.parse_args()
    if args.wake_map and not args.flow_data:
        parser.error("--wake-map requires --flow-data")
    print(json.dumps(compare(args.solver.resolve(), args.run.resolve(), args.pressure_data, args.flow_data, args.wake_map), indent=2))
