"""Decode the public format independently of the C++ writer."""
import json
import math
import pathlib
import struct
import sys

root = pathlib.Path(sys.argv[1])


def decode(entry):
    data = (root / "view" / entry["file"]).read_bytes()
    magic, version, header, nv, ni, nf, r0, r1 = struct.unpack_from("<4s7I", data)
    assert (magic, version, header, r0, r1) == (b"KUTV", 1, 32, 0, 0)
    assert (nv, ni, nf) == (entry["vertices"], entry["indices"], len(entry["fields"]))
    p = list(struct.iter_unpack("<3f", data[32:32 + nv * 12]))
    indices = struct.unpack_from(f"<{ni}I", data, 32 + nv * 12)
    assert ni % 3 == 0 and all(i < nv for i in indices)
    offset = 32 + nv * 12 + ni * 4
    fields = {}
    for f in entry["fields"]:
        assert f["location"] == "vertex" and f["offset_bytes"] == offset and f["unit"]
        n = nv * f["components"]
        values = struct.unpack_from(f"<{n}f", data, offset)
        fields[f["name"]] = values
        if f["components"] == 3:
            values = [math.sqrt(sum(x * x for x in values[i:i + 3])) for i in range(0, n, 3)]
        finite = [v for v in values if math.isfinite(v)]
        if finite:
            assert math.isclose(f["min"], min(finite), rel_tol=1e-6, abs_tol=1e-6)
            assert math.isclose(f["max"], max(finite), rel_tol=1e-6, abs_tol=1e-6)
        else:
            assert f["min"] is None and f["max"] is None
        offset += n * 4
    assert offset == len(data)
    area = 0
    for i in range(0, ni, 3):
        a, b, c = [p[n] for n in indices[i:i + 3]]
        u, v = [b[j] - a[j] for j in range(3)], [c[j] - a[j] for j in range(3)]
        cross = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]]
        normal = entry.get("normal", {"xmin": [-1, 0, 0], "xmax": [1, 0, 0], "ymin": [0, -1, 0], "ymax": [0, 1, 0], "zmin": [0, 0, -1], "zmax": [0, 0, 1]}.get(entry.get("name")))
        assert normal is not None and sum(x * n for x, n in zip(cross, normal)) > 0
        area += math.sqrt(sum(x * x for x in cross)) / 2
    return p, fields, area


manifest = json.loads((root / "view/manifest.json").read_text())
assert manifest["v"] == 1 and manifest["format"]["version"] == 1 and manifest["reference"]
assert manifest["bounds"] == [[0, 0, 0], [2, 1, 1]]
assert len(manifest["zones"]) == 6 and len(manifest["slices"]) == 5
for zone in manifest["zones"]:
    p, fields, area = decode(zone)
    expected_area = 1 if zone["name"] in ("xmin", "xmax") else 2
    assert math.isclose(area, expected_area)
    assert "pressure" in fields and "pressure@delta" in fields
    assert all(v == 2 for v in fields["pressure@delta"])
    for position, pressure in zip(p, fields["pressure"]):
        assert pressure == {0: 10, 1: 15, 2: 20}[position[0]], (zone, position, pressure)
    if zone["type"] == "wall":
        assert {"wall_shear", "wall_shear_magnitude", "y_plus"} <= fields.keys()
        assert "velocity" not in fields
    else:
        assert "velocity" in fields and "wall_shear" not in fields
for n, entry in enumerate(manifest["slices"]):
    p, fields, area = decode(entry)
    assert math.isclose(area, [1, 2, 1, 0, math.sqrt(2)][n], rel_tol=1e-6), (n, area)
    assert all(v == 2 for v in fields["pressure@delta"])
    assert all(v == 1 for v in fields["velocity@delta"])
    if n in (0, 2):
        assert all(v == 10 for v in fields["pressure"])
    if n == 0:
        assert all(v[0] == .5 for v in p)
        assert {v[1:] for v in p} >= {(0, 0), (0, 1), (1, 0), (1, 1)}
print("View format: header, geometry, winding-compatible areas, per-vertex values, ranges, deltas and replacement passed")
