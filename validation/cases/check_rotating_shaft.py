"""Check the packaged Taylor–Couette example against the steady circular flow."""
import argparse
import json
import math
import subprocess
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path


def check(directory):
    root = ET.parse(directory / "result.vtu")
    coords = list(map(float, root.find(".//Points/DataArray").text.split()))
    points = list(zip(coords[::3], coords[1::3], coords[2::3], strict=True))
    cells = list(map(int, root.find(".//Cells/DataArray[@Name='connectivity']").text.split()))
    offsets = list(map(int, root.find(".//Cells/DataArray[@Name='offsets']").text.split()))
    velocity = list(map(float, root.find(".//CellData/DataArray[@Name='u']").text.split()))
    pressure = list(map(float, root.find(".//CellData/DataArray[@Name='p']").text.split()))
    ri, ro, omega, mu, rho = 0.1, 0.2, 10, 0.01, 1.0
    a = -omega * ri**2 / (ro**2 - ri**2)
    b = -a * ro**2
    pressure_samples = []
    square_error = square_reference = 0.0
    start = 0
    for i, end in enumerate(offsets):
        ring = sorted({points[node][:2] for node in cells[start:end]})
        start = end
        if len(ring) != 4:
            raise ValueError("Expected an extruded quadrilateral")
        xm = sum(x for x, y in ring) / 4
        ym = sum(y for x, y in ring) / 4
        ring.sort(key=lambda p: math.atan2(p[1] - ym, p[0] - xm))
        area2 = cx = cy = 0.0
        for (x, y), (xx, yy) in zip(ring, ring[1:] + ring[:1]):
            cross = x * yy - xx * y
            area2 += cross
            cx += (x + xx) * cross
            cy += (y + yy) * cross
        x, y = cx / (3 * area2), cy / (3 * area2)
        r_squared = x*x + y*y
        angular_velocity = a + b / r_squared
        p_exact = rho * (0.5 * a*a * r_squared + a*b * math.log(r_squared / ri**2) - 0.5 * b*b / r_squared)
        pressure_samples.append((area2, p_exact, pressure[i] - p_exact))
        u, v, w = velocity[3*i:3*i+3]
        square_error += area2 * ((u + angular_velocity*y)**2 + (v - angular_velocity*x)**2 + w*w)
        square_reference += area2 * angular_velocity**2 * (x*x + y*y)
    area = sum(w for w, _, _ in pressure_samples)
    p_mean = sum(w * p for w, p, _ in pressure_samples) / area
    gauge_offset = sum(w * delta for w, _, delta in pressure_samples) / area
    p_error = sum(w * (delta - gauge_offset)**2 for w, _, delta in pressure_samples)
    p_variation = sum(w * (p - p_mean)**2 for w, p, _ in pressure_samples)
    zones = json.loads((directory / "loads.json").read_text())["zones"]
    exact = 4 * math.pi * mu * omega * ri**2 * ro**2 / (ro**2 - ri**2)
    errors = {"relative_velocity_l2": math.sqrt(square_error / square_reference),
              "relative_pressure_l2": math.sqrt(p_error / p_variation)}
    for name, target in [("shaft", -exact), ("housing", exact)]:
        wall = next(z for z in zones if z["name"] == name)
        torque = wall["viscous_moment"][2] + wall["pressure_moment"][2]
        errors[f"{name}_torque"] = torque
        errors[f"{name}_relative_error"] = abs(torque / target - 1)
    return errors


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("solver", type=Path)
    parser.add_argument("--backend", default="cpu")
    args = parser.parse_args()
    case = Path(__file__).resolve().parents[2] / "app/resources/examples/rotating-shaft.kutta.json"
    with tempfile.TemporaryDirectory() as temp:
        out = Path(temp)
        subprocess.run([args.solver.resolve(), "solve", case, "--backend", args.backend, "--out", out], check=True)
        summary = json.loads((out / "summary.json").read_text())
        if summary["verdict"] != "converged" or summary["convergence"]["residuals_waived"]:
            raise SystemExit("Rotating-shaft residual criteria missed")
        result = check(out)
        print(json.dumps(result, indent=2))
        if not (result["relative_velocity_l2"] < 0.01 and result["relative_pressure_l2"] < 0.01
                and result["shaft_relative_error"] < 0.02 and result["housing_relative_error"] < 0.02):
            raise SystemExit("Rotating-shaft accuracy target missed")
