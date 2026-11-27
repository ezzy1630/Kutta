"""Match the published Windsor pressure taps to Kutta's yawed wall sections."""
import csv
import math
from pathlib import Path


def sections(experiment, drag, side, walls):
    with (experiment / "SquareBack_Pressure_Tapping_Map.csv").open(newline="") as src:
        taps = list(csv.DictReader(src))
    with (experiment / "SquareBack_NW_2_5yaw_Pressures_Averages.csv").open(newline="") as src:
        header, *records = csv.reader(src)
    mean = next(row for row in records if row[0] == "MeanCp")
    planes = {}
    for tap, label, cp in zip(taps, header[1:], mean[1:], strict=True):
        number = tap["Tapping Number"]
        if label != f"Tapping Number {number}":
            raise ValueError("Tap identities differ between the coordinate map and measurements")
        group = int(tap["Group Number"])
        x, y, z = (float(tap[f"{axis} (m)"]) for axis in "xyz")
        point = (x*drag[0] + y*side[0], x*drag[1] + y*side[1], z)
        key = (group, y if group == 1 else z)
        if key not in planes:
            probe = {"name": f"Tap group {group}, section {key[1]:g} m",
                     "type": "section", "zones": [walls[1] if group == 4 else walls[0]],
                     "origin": list(point), "normal": side if group == 1 else [0, 0, 1],
                     "fields": ["pressure_coefficient"]}
            planes[key] = (probe, [])
        planes[key][1].append((number, group, point, float(cp)))
    return list(planes.values())


def compare_sections(view: Path, planes):
    results = []
    for i, (_, taps) in enumerate(planes, start=1):
        with (view / f"probe-{i}.csv").open(newline="") as src:
            rows = list(csv.DictReader(src))
        segments = []
        for a, b in zip(rows, rows[1:]):
            if b["valid"] != "1" or a["valid"] != "1" or float(b["arc_length"]) <= float(a["arc_length"]):
                continue
            start = tuple(float(a[k]) for k in "xyz")
            delta = tuple(float(b[k]) - s for k, s in zip("xyz", start))
            length2 = sum(d*d for d in delta)
            if length2:
                segments.append((start, delta, length2, float(a["pressure_coefficient"]), float(b["pressure_coefficient"])))
        if not segments:
            raise ValueError(f"Empty section at tap {taps[0][0]}")
        for number, group, point, measured in taps:
            nearest = (math.inf, math.nan)
            for start, delta, length2, ca, cb in segments:
                offset = [p - s for p, s in zip(point, start)]
                t = min(1., max(0., sum(o*d for o, d in zip(offset, delta)) / length2))
                distance = math.sqrt(sum((o - t*d)**2 for o, d in zip(offset, delta)))
                if distance < nearest[0]:
                    nearest = distance, ca + t*(cb - ca)
            distance, predicted = nearest
            results.append({"tap": number, "group": group, "tunnel_position_m": point,
                            "projection_distance_m": distance,
                            "measured_Cp": measured, "predicted_Cp": predicted, "error": predicted - measured})
    by_group = {}
    for group in sorted({p["group"] for p in results}):
        sampled = [p for p in results if p["group"] == group]
        by_group[group] = {"samples": len(sampled),
                           "rms_error_Cp": math.sqrt(sum(p["error"]**2 for p in sampled) / len(sampled)),
                           "max_projection_distance_m": max(p["projection_distance_m"] for p in sampled)}
    return {"groups": by_group, "taps": results}
