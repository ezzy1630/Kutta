"""Real Fluent smoke coverage; analytic acceptance lives in the C++ tests."""
import json
import math
import pathlib
import struct
import subprocess
import sys

cli, fixture, run = sys.argv[1], pathlib.Path(sys.argv[2]), pathlib.Path(sys.argv[3])


def command(*args):
    p = subprocess.run([cli, *map(str, args)], capture_output=True, text=True, timeout=30)
    assert p.returncode == 0, p.stdout + p.stderr
    return p.stdout


command("solve", fixture / "mixing_elbow.cas.h5", "--data", fixture / "mixing_elbow.dat.h5",
        "--approximate", "wfgc", "--backend", "cpu", "--iterations", 1, "--out", run)
command("export-view", run)
base = json.loads((run / "view/manifest.json").read_text())
wall = next(z for z in base["zones"] if z["type"] == "wall")
file = (run / "view" / wall["file"]).read_bytes()
points = list(struct.iter_unpack("<3f", file[32:32 + wall["vertices"] * 12]))
lo = [min(p[a] for p in points) for a in range(3)]
hi = [max(p[a] for p in points) for a in range(3)]
axis = max(range(3), key=lambda a: hi[a] - lo[a])
normal = [int(a == axis) for a in range(3)]
origin = [(lo[a] + hi[a]) / 2 for a in range(3)]
scale = math.sqrt(sum((hi[a] - lo[a]) ** 2 for a in range(3)))
inlet = next(z for z in base["zones"] if z["type"] == "velocity-inlet")
seed_spec = {"v": 1, "name": "inlet", "seeds": {"type": "surface", "zones": [inlet["name"]], "count": 2, "offset": scale * 1e-5},
             "max_length": scale * .05, "max_steps": 40, "fields": ["velocity_magnitude", "pressure"]}
seed_file = run / "seeds.json"
seed_file.write_text(json.dumps(seed_spec))
command("export-view", run, "--streamlines", seed_file)
traced = json.loads((run / "view/manifest.json").read_text())
trace = traced["streamlines"][0]
assert trace["vertices"] > 2 and any(line["length"] > 0 for line in trace["lines"]), trace
positions = (run / "view" / trace["file"]).read_bytes()
seed = struct.unpack_from("<3f", positions, 32)
end = [seed[a] + scale * 1e-6 * normal[a] for a in range(3)]
pressure = next(f for f in base["fields"] if f["name"] == "pressure")
level = (pressure["min"] + pressure["max"]) / 2
post = {"v": 1, "streamlines": [{k: v for k, v in seed_spec.items() if k != "v"}],
        "isosurfaces": [{"name": "pressure-level", "field": "pressure", "value": level},
                        {"name": "q-level", "field": "q_criterion", "value": 0}],
        "probes": [{"name": "line", "type": "line", "start": list(seed), "end": end, "count": 7,
                    "fields": ["pressure", "total_pressure_coefficient", "velocity_magnitude", "k", "omega"]},
                   {"name": "section", "type": "section", "zones": [wall["name"]], "origin": origin, "normal": normal,
                    "fields": ["pressure_coefficient", "skin_friction", "y_plus"]}],
        "integrals": [{"name": "all-zones", "fields": ["pressure", "total_pressure_coefficient"]}]}
post_file = run / "post.json"
post_file.write_text(json.dumps(post))
command("export-view", run, "--post", post_file)
manifest = json.loads((run / "view/manifest.json").read_text())
assert manifest["reference"] and "APPROXIMATE RESULT" in manifest["provenance"]
assert len(manifest["streamlines"]) == 1 and len(manifest["isosurfaces"]) == 2 and len(manifest["probes"]) == 2
assert manifest["isosurfaces"][0]["indices"] > 0
assert manifest["probes"][1]["vertices"] > 0
assert len(manifest["reports"][0]["zones"]) == len(base["zones"])
fields = {f["name"] for f in manifest["fields"]}
assert {"pressure_coefficient", "total_pressure_coefficient", "vorticity", "q_criterion", "helicity", "skin_friction", "velocity_x", "velocity_y", "velocity_z"} <= fields
assert {"q_criterion@fluent", "q_criterion@delta", "total_pressure_coefficient@fluent", "total_pressure_coefficient@delta"} <= fields
for kind in ["streamlines", "isosurfaces", "probes"]:
    for entry in manifest[kind]:
        data = (run / "view" / entry["file"]).read_bytes()
        assert data[:4] in (b"KUTV", b"KUTL")
        for f in entry["fields"]:
            count = entry["vertices"] * f["components"]
            values = struct.unpack_from(f"<{count}f", data, f["offset_bytes"])
            assert not any(math.isinf(x) for x in values)
for report in manifest["reports"][0]["zones"]:
    assert report["area"] > 0 and math.isfinite(report["mass_flow"])
    assert report["force"] is not None and report["moment"] is not None
    assert all(math.isfinite(x) for x in report["force"]["total"] + report["moment"]["total"])
print("Fluent mixing elbow: streamlines, pressure/Q iso-surfaces, line/section probes, all derived/reference/delta fields and all-zone aero reports passed")
