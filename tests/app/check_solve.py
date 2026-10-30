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


def check_outputs(directory, data=False, iterations=21):
    summary = json.loads((directory / "summary.json").read_text())
    assert summary["case_sha256"] == hashlib.sha256((root / "tiny.cas.h5").read_bytes()).hexdigest()
    assert summary["backend"] == "cpu" and summary["device"]
    assert summary["cells"] == 4 and summary["iterations"] == iterations
    assert summary["wall_time_seconds"] > 0 and summary["setup_time_seconds"] > 0
    assert summary["converged"] is False and summary["stop_reason"] == "iteration_limit"
    assert summary["convergence_criteria_available"] is False
    assert summary["peak_device_memory_bytes"] is None
    assert summary["kutta_version"] and summary["kutta_commit"]
    expected = hashlib.sha256((root / "tiny.dat.h5").read_bytes()).hexdigest() if data else None
    assert summary["data_sha256"] == expected
    piece = ET.parse(directory / "result.vtu").find("UnstructuredGrid/Piece")
    assert int(piece.attrib["NumberOfCells"]) == 4
    arrays = {a.attrib["Name"]: list(map(float, a.text.split())) for a in piece.findall("CellData/DataArray")}
    assert set(arrays) == {"u", "p", "k", "omega", "mu_t"}
    for key, values in arrays.items():
        assert len(values) == (12 if key == "u" else 4)
        assert all(math.isfinite(v) for v in values)
    with (directory / "residuals.csv").open(newline="") as f:
        rows = list(csv.DictReader(f))
    assert [int(r["iteration"]) for r in rows] == ([20, 21] if iterations == 21 else [iterations])
    for row in rows:
        assert all(math.isfinite(float(v)) and float(v) >= 0 for v in row.values())
        assert float(row["seconds_per_iteration"]) > 0
    with (directory / "forces.csv").open(newline="") as f:
        rows = list(csv.DictReader(f))
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
    return arrays


case = root / "tiny.cas.h5"
log = run(case, "--backend", "cpu", "--iterations", 21, "--out", root / "standard")
assert all(s in log for s in ["continuity", "x-velocity", "y-velocity", "z-velocity", "omega", "sec/iter"])
check_outputs(root / "standard")
run(case, "--backend", "cpu", "--iterations", 1, "--out", root / "cold")
cold = check_outputs(root / "cold", iterations=1)
run(case, "--data", root / "tiny.dat.h5", "--backend", "cpu", "--iterations", 1, "--out", root / "seeded")
seeded = check_outputs(root / "seeded", data=True, iterations=1)
assert max(map(abs, cold["u"])) > 1e-6 and max(map(abs, seeded["u"])) < 1e-10
assert max(seeded["k"]) > 5 * max(cold["k"]), "--data did not seed turbulence fields"
error = run(root / "unsupported.cas.h5", "--backend", "cpu", "--out", root / "rejected", success=False)
assert "recon/cell-lsf?" in error and "curvature-correction?" in error, error
assert not (root / "rejected").exists()
for args in [("--iterations", "0"), ("--iterations", "1.5"), ("--iterations", "999999999999999"), ("--backend", "invalid"), ("--out",), ("--iterations", "1", "--iterations", "2")]:
    run(case, *args, success=False)
print("Solve CLI: generated CFF/data, residuals, VTU, forces, SHA256, and rejection checks passed")
