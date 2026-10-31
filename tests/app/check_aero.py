"""Real Fluent case: convergence stop, loads.json against forces.csv and kutta aero."""
import csv
import json
import math
import pathlib
import subprocess
import sys

cli, fixture, root = sys.argv[1], pathlib.Path(sys.argv[2]), pathlib.Path(sys.argv[3])


def command(*args, code=0):
    p = subprocess.run([cli, *map(str, args)], capture_output=True, text=True, timeout=120)
    assert p.returncode == code, p.stdout + p.stderr
    return p.stdout


# Only continuity is checked, with a criterion every scaled residual meets, so
# the first monitor sample converges: a deterministic stop before the limit.
overrides = {"run/monitor-interval": 2, "monitors/check-convergence": True, "monitors/residual/continuity/absolute-criterion": 10,
             **{f"monitors/residual/{eq}/check": False for eq in ("x-velocity", "y-velocity", "z-velocity", "k", "omega")}}
setup = root / "setup.json"
setup.write_text(json.dumps({"v": 1, "overrides": overrides}))
run = root / "run"
events = [json.loads(line) for line in command("solve", fixture / "mixing_elbow.cas.h5", "--approximate", "wfgc", "--backend", "cpu",
                                               "--iterations", 40, "--setup", setup, "--out", run, "--progress", "jsonl").splitlines()]
assert events[0]["event"] == "setup" and events[0]["modified"] is False and events[0]["run_controls"] == overrides
assert events[0]["convergence"] == {"enabled": True, "condition": "all", "residuals": {"continuity": 10}, "reports": []}
# Iteration 1 is reported (already monitored) but is not a convergence sample.
early, sample = [e for e in events if e["event"] == "iteration"]
assert early["iter"] == 1 and "convergence" not in early and early["reports"] == {}, early
assert sample["iter"] == 2 and sample["convergence"]["met"], sample
[condition] = sample["convergence"]["conditions"]
assert condition["kind"] == "residuals" and list(condition["equations"]) == ["continuity"]
assert condition["equations"]["continuity"]["value"] == sample["residuals"]["continuity"] <= 10
done = events[-1]
assert done["event"] == "done" and done["converged"] is True and done["stop_reason"] == "converged" and done["iterations"] == 2
assert done["convergence"]["met_by"] == ["residuals"] and done["convergence"]["iteration"] == 2
summary = json.loads((run / "summary.json").read_text())
assert summary["converged"] is True and summary["stop_reason"] == "converged" and summary["convergence"] == done["convergence"]
assert summary["convergence_criteria_available"] is True and summary["iterations"] == 2

# Per-zone loads equal the solver's wall force rows up to float summation.
loads = json.loads((run / "loads.json").read_text())
assert loads["iterations"] == 2 and loads["approximate"] is True and "APPROXIMATE RESULT" in loads["provenance"]
assert loads["reference"] == {"area": 1, "density": 1.225, "velocity": 1, "length": 1, "pressure": 0}, loads["reference"]
assert loads["report_zones"] == [] and loads["units"] == {"force": "N", "moment": "N m", "area": "m^2"}
with (run / "forces.csv").open(newline="") as f:
    walls = [r for r in csv.DictReader(line for line in f if not line.startswith("#")) if r["kind"] == "wall"]
assert [z["name"] for z in loads["zones"]] == [r["name"] for r in walls] and len(walls) >= 2
worst = 0
for zone, row in zip(loads["zones"], walls):
    assert zone["zone_id"] == int(row["zone_id"]) and zone["area"] > 0
    total = [zone["pressure"][i] + zone["viscous"][i] for i in range(3)]
    expected = [float(row[f"total_{axis}"]) for axis in "xyz"]
    worst = max(worst, math.dist(total, expected) / math.hypot(*expected))
assert worst < 1e-5, worst

# Aero coefficients from loads.json only; the default car is every wall zone.
[line] = command("aero", run, "--json").splitlines()
event = json.loads(line)
assert event["event"] == "aero" and event["v"] == 1 and event["iterations"] == 2 and event["approximate"] is True
stored = json.loads((run / "aero.json").read_text())
assert stored == {k: v for k, v in event.items() if k != "event"}
assert [r["name"] for r in event["rows"]] == [z["name"] for z in loads["zones"]] + ["Total"]
total = event["rows"][-1]
assert total["kind"] == "total" and total["zones"] == [z["name"] for z in loads["zones"]]
drag = sum(z["pressure"][0] + z["viscous"][0] for z in loads["zones"])
assert math.isclose(total["force"]["total"]["drag"], drag, rel_tol=1e-12)
assert math.isclose(total["coefficients"]["cd"], drag / (.5 * 1.225), rel_tol=1e-12)
assert event["balance"]["status"] == "needs_vehicle" and event["reference"]["dynamic_pressure"] == .5 * 1.225
assert not [p for p in run.iterdir() if p.name.startswith(".aero")], "staging file leaked"

first, second = (z["name"] for z in loads["zones"][:2])
spec = root / "aero.json"
spec.write_text(json.dumps({"v": 1, "groups": [{"name": "Pipe", "zones": [first, second]}],
                            "vehicle": {"front_axle": [0, 0, 0], "wheelbase": 1}}))
grouped = json.loads(command("aero", run, "--aero", spec, "--json"))
pipe = next(r for r in grouped["rows"] if r["kind"] == "group")
assert pipe["name"] == "Pipe" and math.isclose(pipe["downforce"], sum(r["downforce"] for r in grouped["rows"][:2]), rel_tol=1e-12, abs_tol=1e-15)
assert grouped["balance"]["status"] == ("defined" if total["downforce"] > 0 else "no_downforce")
human = command("aero", run)
assert "Total" in human and "Wrote" in human

# Invalid specs are usage errors; a run without loads.json is a failure.
spec.write_text(json.dumps({"v": 1, "car_zones": ["no-such-zone"]}))
[error] = map(json.loads, command("aero", run, "--aero", spec, "--json", code=2).splitlines())
assert error["event"] == "error" and "no-such-zone" in error["message"]
spec.write_text('{"v":1,"axes":{"drag":[1,0,0],"lift":[1,0,0]}}')
command("aero", run, "--aero", spec, code=2)
assert json.loads((run / "aero.json").read_text()) == stored, "a rejected spec replaced aero.json"
(root / "old").mkdir()
[error] = map(json.loads, command("aero", root / "old", "--json", code=1).splitlines())
assert error["event"] == "error" and "loads.json" in error["message"]
command("aero", code=2)
command("aero", run, "--bogus", code=2)
print(f"Aero: converged at iteration 2, loads.json vs forces.csv max relative difference {worst:.3g}, aero rows and errors passed")
