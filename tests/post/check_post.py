"""Decode KUTL/KUTV through the actual CLI, without using C++ writer code."""
import copy
import csv
import json
import math
import pathlib
import struct
import subprocess
import sys

cli, root = sys.argv[1], pathlib.Path(sys.argv[2])


def run(spec, *, code=0, flag="--post"):
    path = root / "post.json"
    path.write_text(json.dumps(spec))
    result = subprocess.run([cli, "export-view", str(root), flag, str(path)], capture_output=True, text=True, timeout=20)
    assert result.returncode == code, result.stdout + result.stderr
    assert not list(root.glob(".view-*")), "staging directory leaked"
    if code == 0:
        return json.loads((root / "view/manifest.json").read_text())


def decode(entry):
    data = (root / "view" / entry["file"]).read_bytes()
    magic, version, size, nv, ni, nf, lines, reserved = struct.unpack_from("<4s7I", data)
    assert magic in (b"KUTL", b"KUTV") and version == 1 and size == 32 and reserved == 0
    assert (nv, ni, nf) == (entry["vertices"], entry["indices"], len(entry["fields"]))
    positions = list(struct.iter_unpack("<3f", data[32:32 + nv * 12]))
    assert all(math.isfinite(x) for p in positions for x in p)
    indices = struct.unpack_from(f"<{ni}I", data, 32 + nv * 12)
    assert all(n < nv for n in indices)
    offset = 32 + nv * 12 + ni * 4
    if magic == b"KUTL":
        assert entry["primitive"] == "polylines" and entry["line_count"] == lines
        offsets = struct.unpack_from(f"<{lines + 1}I", data, offset)
        assert offsets[0] == 0 and offsets[-1] == ni and list(offsets) == sorted(offsets)
        offset += (lines + 1) * 4
    else:
        assert entry["primitive"] == "triangles" and ni % 3 == 0 and lines == 0
        offsets = []
    fields = {}
    for f in entry["fields"]:
        assert f["offset_bytes"] == offset and f["location"] == "vertex" and f["unit"]
        count = nv * f["components"]
        values = struct.unpack_from(f"<{count}f", data, offset)
        fields[f["name"]] = values
        if f["components"] == 3:
            values = [math.sqrt(sum(x * x for x in values[i:i + 3])) for i in range(0, count, 3)]
        finite = [v for v in values if math.isfinite(v)]
        if finite:
            assert math.isclose(min(finite), f["min"], rel_tol=2e-6, abs_tol=1e-6)
            assert math.isclose(max(finite), f["max"], rel_tol=2e-6, abs_tol=1e-6)
        else:
            assert f["min"] is None and f["max"] is None
        offset += count * 4
    assert offset == len(data)
    return positions, indices, offsets, fields


trace = {"name": "rake", "seeds": {"type": "rake", "start": [.1, .2, .5], "end": [.1, .8, .5], "count": 3},
         "max_length": 2, "max_steps": 500, "fields": ["pressure", "velocity", "q_criterion"]}
spec = {"v": 1, "streamlines": [trace],
        "isosurfaces": [{"name": "iso", "field": "phi", "value": .43, "fields": ["pressure"]}],
        "probes": [{"name": "line", "type": "line", "start": [0, .5, .5], "end": [1, .5, .5], "count": 11, "fields": ["phi", "pressure"]},
                   {"name": "section", "type": "section", "zones": ["xmax"], "origin": [0, .5, 0], "normal": [0, 1, 0], "fields": ["pressure_coefficient"]}],
        "integrals": [{"zones": ["xmax"], "fields": ["pressure", "skin_friction"]}]}
manifest = run(spec)
assert manifest["v"] == 1
for key in ["streamlines", "isosurfaces", "probes", "reports", "available_fields", "reference_values"]:
    assert key in manifest
for key in ["streamlines", "isosurfaces", "probes"]:
    for entry in manifest[key]:
        decode(entry)
entry = manifest["streamlines"][0]
p, indices, offsets, fields = decode(entry)
assert len(entry["lines"]) == 3 and all(x["stop_reason"] == "boundary" for x in entry["lines"])
for a, b in zip(offsets, offsets[1:]):
    assert math.isclose(p[indices[a]][0], .1, abs_tol=1e-6)
    assert math.isclose(p[indices[b - 1]][0], 1, abs_tol=1e-6)
p, _, _, _ = decode(manifest["isosurfaces"][0])
assert all(math.isclose(x[0], .43, abs_tol=1e-6) for x in p)
with (root / "view" / manifest["probes"][0]["csv"]).open() as f:
    rows = list(csv.DictReader(f))
assert len(rows) == 11 and math.isclose(float(rows[4]["phi"]), .4, abs_tol=1e-6)
reports = json.loads((root / "view" / manifest["reports_file"]).read_text())
assert reports["v"] == 1 and reports["reports"] == manifest["reports"]
assert reports["reports"][0]["total"]["mass_flow_source"] == "solver_face_flux"
available = {f["name"] for f in manifest["available_fields"]}
assert {"total_pressure_coefficient", "vorticity", "q_criterion", "helicity", "skin_friction", "velocity_x", "velocity_y", "velocity_z"} <= available

previous = (root / "view/manifest.json").read_bytes()
for bad, code in [({"v": 2}, 2), ({"v": 1, "unknown": 1}, 2),
                  ({"v": 1, "probes": [{"start": [0, 0, 0], "end": [1, 0, 0], "count": 0}]}, 2),
                  ({"v": 1, "isosurfaces": [{"field": "unknown", "value": 0}]}, 1),
                  ({"v": 1, "integrals": [{"zones": ["not-a-zone"]}]}, 1),
                  ({"v": 1, "reference": {"density": 3}}, 1)]:
    run(bad, code=code)
    assert (root / "view/manifest.json").read_bytes() == previous

# The requested shortcut, and true time interpolation through files, not test
# callback velocities. Displacement over t=[0,.1] is .1 + .5*.1^2.
pathline = copy.deepcopy(trace)
pathline.update({"v": 1, "mode": "pathline", "name": "unsteady", "max_time": .1,
                 "seeds": {"type": "rake", "start": [.1, .5, .5], "end": [.1, .5, .5], "count": 1},
                 "frames": [{"time": 0, "run": "."}, {"time": 2, "run": "frame1"}]})
manifest = run(pathline, flag="--streamlines")
assert not manifest["isosurfaces"] and not list((root / "view").glob("iso-*.bin"))
p, _, _, fields = decode(manifest["streamlines"][0])
assert math.isclose(p[-1][0], .205, abs_tol=1e-6), p[-1]
assert math.isclose(fields["time"][-1], .1, abs_tol=1e-6)
print("Post CLI: KUTV/KUTL layout, fields/ranges, pathline frame interpolation, probes/CSV, reports and atomic failures passed")
