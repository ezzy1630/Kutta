"""Real CLI checks; uses only Python's standard library on every CI platform."""
import csv
import hashlib
import json
import math
import pathlib
import subprocess
import sys
import xml.etree.ElementTree as ET

cli, root = sys.argv[1], pathlib.Path(sys.argv[2])


def run(*args, success=True):
    result = subprocess.run([cli, "solve", *map(str, args)], capture_output=True, text=True)
    assert (result.returncode == 0) == success, result.stdout + result.stderr
    return result.stdout + result.stderr


def csv_rows(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(line for line in f if not line.startswith("#")))


def check_outputs(directory, data=False, iterations=21, case_name="tiny.cas.h5", approximate=False):
    summary = json.loads((directory / "summary.json").read_text())
    assert summary["case_sha256"] == hashlib.sha256((root / case_name).read_bytes()).hexdigest()
    assert summary["backend"] == "cpu" and summary["device"]
    assert summary["cells"] == 4 and summary["iterations"] == iterations
    assert summary["wall_time_seconds"] > 0 and summary["setup_time_seconds"] > 0
    assert summary["converged"] is False and summary["stop_reason"] == "iteration_limit"
    assert summary["convergence_criteria_available"] is True
    assert summary["convergence"] == {"enabled": False, "condition": "all", "residuals": {}, "reports": [], "met_by": [], "iteration": None}
    assert summary["peak_device_memory_bytes"] is None
    assert summary["kutta_version"] and summary["kutta_commit"]
    assert summary["approximate"] is approximate
    if approximate:
        [info] = summary["approximations"]
        assert info["name"] == "wfgc" and info["setting"] == "recon/cell-lsf?"
        assert info["fluent_value"] == "#t"
        assert "WFGC off, #f" in info["kutta_action"]
        assert info["effect_documentation"] == "src/kutta/fv/VERIFICATION.md#v6-initialization-comparison"
        assert "3.366748%" in info["measured_effect"] and "not isolated WFGC error" in info["measured_effect"]
        assert info["short"] == "WFGC off"
        notice = summary["approximation_notice"]
        assert notice.startswith("APPROXIMATE RESULT")
        for key in ("setting", "fluent_value", "kutta_action", "effect_documentation", "measured_effect"):
            assert info[key] in notice
        for name in ("residuals.csv", "forces.csv"):
            assert (directory / name).read_text().splitlines()[0] == "# " + notice
        # Parse the header as a comment, proving that the warning is attached to
        # the VTU itself, outside the cell fields and independent of summary.json.
        comments = [e.text for _, e in ET.iterparse(directory / "result.vtu", events=("comment",))]
        assert any("APPROXIMATE RESULT" in c and info["setting"] in c and info["effect_documentation"] in c for c in comments)
    else:
        assert summary["approximations"] == [] and summary["approximation_notice"] == ""
        assert not any("APPROXIMATE RESULT" in (directory / name).read_text() for name in ("residuals.csv", "forces.csv", "result.vtu"))
    expected = hashlib.sha256((root / "tiny.dat.h5").read_bytes()).hexdigest() if data else None
    assert summary["data_sha256"] == expected
    piece = ET.parse(directory / "result.vtu").find("UnstructuredGrid/Piece")
    assert int(piece.attrib["NumberOfCells"]) == 4
    arrays = {a.attrib["Name"]: list(map(float, a.text.split())) for a in piece.findall("CellData/DataArray")}
    assert set(arrays) == {"u", "p", "k", "omega", "mu_t"}
    for key, values in arrays.items():
        assert len(values) == (12 if key == "u" else 4)
        assert all(math.isfinite(v) for v in values)
    rows = csv_rows(directory / "residuals.csv")
    # Iterations 1-5 (already monitored for Fluent's continuity scale) and the samples.
    assert [int(r["iteration"]) for r in rows] == ([1, 2, 3, 4, 5, 20, 21] if iterations == 21 else [iterations])
    for row in rows:
        assert all(math.isfinite(float(v)) and float(v) >= 0 for v in row.values())
        assert float(row["seconds_per_iteration"]) > 0
    rows = csv_rows(directory / "forces.csv")
    assert len(rows) == 4
    assert {(r["kind"], r["name"]) for r in rows} == {("wall", "shell"), ("report", "drag"), ("report", "lift"), ("report", "force")}
    for row in rows:
        assert None not in row, row
        for axis in "xyz":
            assert math.isclose(float(row[f"pressure_{axis}"]) + float(row[f"viscous_{axis}"]), float(row[f"total_{axis}"]), abs_tol=1e-12)
        if row["name"] == "drag":
            assert math.isclose(float(row["coefficient"]), float(row["total_projection"]) / 6, abs_tol=1e-12)
        else:
            assert row["coefficient"] == ""
    # loads.json integrates the same wall data per zone as the forces.csv rows.
    loads = json.loads((directory / "loads.json").read_text())
    assert loads["v"] == 1 and loads["origin"] == [0, 0, 0] and loads["iterations"] == iterations
    assert loads["reference"] == {"area": 3, "length": 1, "density": 1, "velocity": 2, "pressure": 0}
    assert loads["report_zones"] == ["shell"] and loads["approximate"] is approximate and loads["modified"] is False
    assert bool(loads["provenance"]) is approximate
    [wall] = [r for r in rows if r["kind"] == "wall"]
    [zone] = loads["zones"]
    assert zone["name"] == "shell" and zone["zone_id"] == int(wall["zone_id"]) and math.isclose(zone["area"], 10)
    # forces.csv sums float face products on the device; compare at the scale
    # of the largest component, since a closed box nearly cancels some.
    scale = max(abs(float(wall[f"{part}_{axis}"])) for part in ("pressure", "viscous") for axis in "xyz")
    for part in ("pressure", "viscous"):
        for i, axis in enumerate("xyz"):
            assert math.isclose(zone[part][i], float(wall[f"{part}_{axis}"]), abs_tol=1e-5 * scale), (part, zone, wall)
    return arrays


wfgc = root / "wfgc.cas.h5"
error = run(wfgc, "--backend", "cpu", "--out", root / "strict-wfgc", success=False)
assert "--approximate wfgc" in error, error
assert not (root / "strict-wfgc").exists()
log = run(wfgc, "--approximate", "wfgc", "--backend", "cpu", "--iterations", 1, "--out", root / "approximate")
assert log.count("APPROXIMATE RESULT") == 2, log
assert log.startswith("APPROXIMATE RESULT") and log.rstrip().splitlines()[-1].startswith("APPROXIMATE RESULT")
approximate_fields = check_outputs(root / "approximate", iterations=1, case_name="wfgc.cas.h5", approximate=True)

case = root / "tiny.cas.h5"
log = run(case, "--backend", "cpu", "--iterations", 21, "--out", root / "standard")
assert all(s in log for s in ["continuity", "x-velocity", "y-velocity", "z-velocity", "omega", "sec/iter"])
assert "Convergence checking is off (monitors/check-convergence)" in log and "Reached the iteration limit (21 iterations)." in log, log
check_outputs(root / "standard")
run(case, "--approximate", "wfgc", "--backend", "cpu", "--iterations", 1, "--out", root / "unused-opt-in")
check_outputs(root / "unused-opt-in", iterations=1)
run(case, "--backend", "cpu", "--iterations", 1, "--out", root / "cold")
cold = check_outputs(root / "cold", iterations=1)
assert approximate_fields == cold, "WFGC opt-in changed more than the WFGC-on setting"
run(case, "--data", root / "tiny.dat.h5", "--backend", "cpu", "--iterations", 1, "--out", root / "seeded")
seeded = check_outputs(root / "seeded", data=True, iterations=1)
assert max(map(abs, cold["u"])) > 1e-6 and max(map(abs, seeded["u"])) < 1e-10
assert max(seeded["k"]) > 5 * max(cold["k"]), "--data did not seed turbulence fields"
log = run(root / "hybrid.cas.h5", "--backend", "cpu", "--iterations", 1, "--out", root / "hybrid")
assert "Initialization: hybrid" in log
hybrid = check_outputs(root / "hybrid", iterations=1, case_name="hybrid.cas.h5")
assert max(map(abs, hybrid["u"])) < 1e-10, "hybrid dispatch used standard velocity defaults"
run(root / "no-reports.cas.h5", "--backend", "cpu", "--iterations", 1, "--out", root / "no-reports")
with (root / "no-reports" / "forces.csv").open(newline="") as f:
    rows = list(csv.DictReader(f))
assert len(rows) == 1 and rows[0]["kind"] == "wall"
error = run(root / "unsupported.cas.h5", "--backend", "cpu", "--out", root / "rejected", success=False)
assert all(key in error for key in ["recon/cell-lsf?", "curvature-correction?", "pseudo-time-step", "omega/default", "sst-pk-factor"]), error
assert not (root / "rejected").exists()
error = run(root / "unsupported.cas.h5", "--approximate", "wfgc", "--backend", "cpu", "--out", root / "still-rejected", success=False)
assert all(key in error for key in ["curvature-correction?", "pseudo-time-step", "omega/default", "sst-pk-factor"]), error
assert "recon/cell-lsf?" not in error, error
assert not (root / "still-rejected").exists()
for value in ["all", "curvature-correction?", "wfgc,curvature-correction?", "", ",wfgc", "wfgc,", "wfgc,wfgc", "WFGC"]:
    error = run(wfgc, "--approximate", value, "--out", root / "invalid-approximation", success=False)
    assert "--approximate" in error
    assert not (root / "invalid-approximation").exists()
for args in [("--iterations", "0"), ("--iterations", "1.5"), ("--iterations", "999999999999999"), ("--backend", "invalid"), ("--out",), ("--iterations", "1", "--iterations", "2"), ("--approximate", "wfgc", "--approximate", "wfgc")]:
    run(case, *args, success=False)
print("Solve CLI: generated CFF/data, residuals, VTU, forces, SHA256, and rejection checks passed")
