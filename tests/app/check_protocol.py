"""Exercise the sidecar boundary; keep stdin open until explicitly stopping."""
import hashlib
import concurrent.futures
import json
import math
import os
import pathlib
import signal
import subprocess
import sys

cli, root = sys.argv[1], pathlib.Path(sys.argv[2])


def command(*args, code=0):
    result = subprocess.run([cli, *map(str, args)], capture_output=True, text=True, timeout=15)
    assert result.returncode == code, result.stdout + result.stderr
    return result


inventory = json.loads(command("devices", "--json").stdout)
assert inventory["v"] == 1
assert sum(d["default"] for d in inventory["devices"]) == 1
assert any(d["backend"] == "cpu" for d in inventory["devices"])
for d in inventory["devices"]:
    assert d["backend"] in ("cpu", "metal", "cuda") and d["name"]
    assert isinstance(d["memory_bytes"], int) and d["memory_bytes"] >= 0
command("devices", "--bogus", code=2)

for name, status in [("tiny", "supported"), ("wfgc", "approximate"), ("unsupported", "unsupported")]:
    case = root / (name + ".cas.h5")
    summary = json.loads(command("inspect", case, "--json").stdout)
    assert summary["v"] == 1 and summary["support"]["status"] == status, summary
    assert summary["file"]["sha256"] == hashlib.sha256(case.read_bytes()).hexdigest()
    assert summary["mesh"] == {"cells": 4, "faces": 20, "nodes": 18,
                               "cell_types": {"tet": 0, "hex": 4, "poly": 0},
                               "bounds": [[0, 0, 0], [2, 1, 1]]}, summary["mesh"]
    assert summary["models"] == {"viscous": "k-omega-sst", "energy": False, "energy_coupled": False}
    assert next(z for z in summary["zones"] if z["name"] == "shell")["faces"] == 16
    assert {r["name"] for r in summary["reports"]} == {"drag", "lift", "force"}
    assert summary["materials"] == [{"name": "water", "density": 1, "viscosity": .01}]
    for issue in summary["support"]["issues"]:
        assert issue["class"] in ("unsupported", "approximable", "informational")
        assert issue["setting"] and issue["label"] and len(issue["explanation"]) > 80
    if name == "wfgc":
        [issue] = [i for i in summary["support"]["issues"] if i["class"] == "approximable"]
        assert issue["approximation"] == "wfgc" and "3.601459%" in issue["explanation"]
        assert "not isolated WFGC error" in issue["explanation"] and issue["evidence"]
    if name == "unsupported":
        assert {"curvature-correction?", "pseudo-time-step", "omega/default", "sst-pk-factor"} <= {i["setting"] for i in summary["support"]["issues"]}
seed = json.loads(command("inspect", root / "tiny.cas.h5", "--data", root / "tiny.dat.h5", "--json").stdout)
assert seed["data"] == {"present": True, "fields": ["SV_K", "SV_O", "SV_P", "SV_U", "SV_V", "SV_W", "SV_WALL_SHEAR", "SV_WALL_YPLUS"], "iterations_hint": 10}
command("inspect", root / "tiny.cas.h5", "--json", "--json", code=2)


def start(name="tiny", *, iterations=21, data=False, approximate=False, suffix="", stdin=subprocess.PIPE):
    out = root / ("jsonl-" + name + suffix)
    args = [cli, "solve", str(root / (name + ".cas.h5")), "--backend", "cpu", "--iterations", str(iterations),
            "--out", str(out), "--progress", "jsonl"]
    if data:
        args += ["--data", str(root / "tiny.dat.h5")]
    if approximate:
        args += ["--approximate", "wfgc"]
    return subprocess.Popen(args, stdin=stdin, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True), out


def finish(p):
    # Read both streams while retaining the stdin writer. communicate() would
    # close stdin and accidentally turn every ordinary run into a stop test.
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        stdout = pool.submit(p.stdout.read)
        stderr = pool.submit(p.stderr.read)
        try:
            code = p.wait(timeout=15)
        except subprocess.TimeoutExpired:
            p.kill()
            p.wait()
            raise
        text, error = stdout.result(), stderr.result()
    assert code == 0, error + text
    events = list(map(json.loads, text.splitlines()))
    assert all(e["v"] == 1 for e in events)
    return events


for name, data, approximate in [("tiny", False, False), ("tiny", True, False), ("hybrid", False, False), ("wfgc", False, True), ("downforce", False, False)]:
    p, out = start(name, data=data, approximate=approximate, suffix="-data" if data else "")
    events = finish(p)
    assert events[0]["event"] == "setup" and events[-1]["event"] == "done", events
    assert events[0]["init"] == ("fluent-data" if data else "hybrid" if name == "hybrid" else "standard")
    assert events[0]["approximations"] == (["wfgc"] if approximate else [])
    iterations = [e for e in events if e["event"] == "iteration"]
    assert [e["iter"] for e in iterations] == [20, 21]
    for e in iterations:
        assert set(e["residuals"]) == {"continuity", "x-velocity", "y-velocity", "z-velocity", "k", "omega"}
        assert all(math.isfinite(v) and v >= 0 for v in e["residuals"].values())
        assert e["sec_per_iter"] > 0 and e["elapsed"] > 0
        assert set(e["reports"]) == {"drag", "downforce" if name == "downforce" else "lift", "force"}
        for r in e["reports"].values():
            assert math.isclose(r["coefficient"], r["force"] / 6, abs_tol=1e-12)
    assert events[-1]["iterations"] == 21 and events[-1]["stop_reason"] == "iteration_limit"
    assert (out / "view-input.bin").is_file()
    command("export-view", out, "--slice", "0.5,0,0,1,0,0")
    manifest = json.loads((out / "view/manifest.json").read_text())
    assert manifest["reference"] is data
    fields = {f["name"] for f in manifest["fields"]}
    assert {"pressure", "pressure_coefficient", "wall_shear", "wall_shear_magnitude", "y_plus", "velocity", "velocity_magnitude", "k", "omega", "viscosity_ratio"} <= fields
    if data:
        assert {"pressure@fluent", "pressure@delta", "velocity@fluent", "velocity@delta", "k@fluent", "k@delta"} <= fields
        assert "viscosity_ratio@fluent" not in fields  # unavailable saved field stays absent
        assert {"wall_shear@fluent", "wall_shear@delta", "wall_shear_magnitude@fluent", "y_plus@fluent"} <= fields
        for zone in manifest["zones"]:
            refs = {f["name"]: f for f in zone["fields"]}
            assert refs["wall_shear_magnitude@fluent"]["min"] == refs["wall_shear_magnitude@fluent"]["max"] == 5
            assert refs["y_plus@fluent"]["min"] == refs["y_plus@fluent"]["max"] == 7
            assert refs["pressure@fluent"]["min"] == refs["pressure@fluent"]["max"] == 2
    assert bool(manifest["provenance"]) is approximate
    p.stdin.close()

# Close stdin and send SIGINT only after setup was observed. Output files and a
# terminal done event must survive both real OS paths. No timing sleeps.
for stop in (["eof", "sigint"] if os.name != "nt" else ["eof"]):
    p, out = start(iterations=1000000, suffix="-stop-" + stop)
    setup = json.loads(p.stdout.readline())
    assert setup["event"] == "setup"
    prefix = []
    if stop == "eof":
        p.stdin.close()
    else:
        # Prove stopping an active iteration stream, rather than just cancelling
        # during initialization before any outer iteration has happened.
        while not prefix or prefix[-1]["event"] != "iteration":
            prefix.append(json.loads(p.stdout.readline()))
        p.send_signal(signal.SIGINT)
    events = prefix + finish(p)
    assert events[-1]["event"] == "done" and events[-1]["stop_reason"] == "stopped", events
    assert 0 <= events[-1]["iterations"] < 1000000
    if stop == "sigint":
        assert events[-1]["iterations"] >= 20
        assert [e for e in events if e["event"] == "iteration"][-1]["iter"] == events[-1]["iterations"]
    summary = json.loads((out / "summary.json").read_text())
    assert summary["stop_reason"] == "stopped" and summary["iterations"] == events[-1]["iterations"]
    for file in ["result.vtu", "forces.csv", "residuals.csv", "view-input.bin"]:
        assert (out / file).is_file(), file
    command("export-view", out)
    if stop != "eof":
        p.stdin.close()

# EOF on the stdin pipe before setup yields a valid initialized, zero-iteration run.
p, out = start(suffix="-early-eof")
p.stdin.close()
events = finish(p)
assert events[-1]["iterations"] == 0 and events[-1]["stop_reason"] == "stopped"
# Stdin that is not a pipe (scripts, CI, nohup) never means stop.
p, out = start(suffix="-devnull", stdin=subprocess.DEVNULL)
events = finish(p)
assert events[-1]["iterations"] == 21 and events[-1]["stop_reason"] == "iteration_limit", events
for name, code, args in [("unsupported", 1, []), ("wfgc", 1, []), ("tiny", 2, ["--iterations", "0"])]:
    result = command("solve", root / (name + ".cas.h5"), "--progress", "jsonl", *args, code=code)
    events = list(map(json.loads, result.stdout.splitlines()))
    assert [e["event"] for e in events] == ["error"]
    assert events[0]["message"] and "issues" in events[0]
    if code == 1:
        assert events[0]["issues"]
command("export-view", out, "--slice", "0,0,0,0,0,0", code=2)
command("export-view", out, "--slice", "0,0,0,1,0", code=2)
print("App protocol: devices, inspect, JSONL event order/reports, EOF and SIGINT passed")

# Overrides are applied before setup and cannot claim editability without a
# bound setter. Exercise the real mass-flux result, not just the echoed JSON.
import csv


def read_fluxes(out):
    with (out / "fluxes.csv").open() as f:
        return {r["name"]: float(r["mass_flow_outward"]) for r in csv.DictReader(line for line in f if not line.startswith("#"))}


def solve_setup(overrides, name, *, reject=False):
    setup_file = root / ("setup-" + name + ".json")
    setup_file.write_text(json.dumps({"v": 1, "overrides": overrides}))
    out = root / ("overridden-" + name)
    p = subprocess.Popen([cli, "solve", str(root / "inlet.cas.h5"), "--backend", "cpu", "--iterations", "1", "--setup", str(setup_file), "--out", str(out), "--progress", "jsonl"],
                         stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if reject:
        stdout, stderr = p.communicate(timeout=15)
        assert p.returncode == 1, stdout + stderr
        events = list(map(json.loads, stdout.splitlines()))
        assert [e["event"] for e in events] == ["error"], events
        assert set(overrides) <= {i["path"] for i in events[0]["issues"]}, events
        assert not out.exists()
        return
    events = finish(p)
    p.stdin.close()
    assert events[0]["overrides"] == overrides and events[0]["modified"] is bool(overrides)
    summary = json.loads((out / "summary.json").read_text())
    assert summary["overrides"] == overrides and summary["modified"] is bool(overrides)
    assert summary["approximate"] is False and summary["approximation_notice"] == ""
    if overrides:
        assert "MODIFIED RESULT" in (out / "result.vtu").read_text()
        assert "MODIFIED RESULT" in (out / "fluxes.csv").read_text()
        command("export-view", out)
        assert "MODIFIED RESULT" in json.loads((out / "view/manifest.json").read_text())["provenance"]
    return out, events


inlet_summary = json.loads(command("inspect", root / "inlet.cas.h5", "--json").stdout)
settings = {s["path"]: s for s in inlet_summary["settings"]}
speed = settings["boundary/inlet/velocity-magnitude"]
assert speed["editable"] and speed["kind"] == "number" and speed["value"] == 1 and speed["unit"] == "m/s"
assert not settings["general/operating-pressure"]["editable"]
assert not settings["initialization/method"]["editable"]
assert settings["run/iterations"]["editable"] and settings["methods/momentum-scheme"]["options"] == ["first-order-upwind", "second-order-upwind"]
source_hash = hashlib.sha256((root / "inlet.cas.h5").read_bytes()).hexdigest()
base, _ = solve_setup({}, "baseline")
changed, events = solve_setup({"boundary/inlet/velocity-magnitude": 2, "run/iterations": 3, "run/monitor-interval": 1}, "inlet")
assert events[-1]["iterations"] == 3
assert [e["iter"] for e in events if e["event"] == "iteration"] == [1, 2, 3]
base_flux, changed_flux = read_fluxes(base), read_fluxes(changed)
assert math.isclose(base_flux["inlet"], -1, abs_tol=1e-6), base_flux
assert math.isclose(changed_flux["inlet"], -2, abs_tol=1e-6), changed_flux
assert not math.isclose(base_flux["outlet"], changed_flux["outlet"], abs_tol=1e-6)
assert hashlib.sha256((root / "inlet.cas.h5").read_bytes()).hexdigest() == source_hash
solve_setup({"missing/setting": 1, "general/operating-pressure": 10, "run/iterations": 1.5,
             "boundary/inlet/velocity-magnitude": -1, "methods/momentum-scheme": "central",
             "run/monitor-interval": "ten", "boundary/inlet/direction": [0, 0, 0]}, "invalid", reject=True)
solve_setup({"run/iterations": True}, "bool-is-not-integer", reject=True)
solve_setup({"boundary/inlet/direction": [1, 0]}, "bad-vector", reject=True)
print(f"Setup overrides: actual inlet mass flux {base_flux['inlet']} -> {changed_flux['inlet']} kg/s; rejections and provenance passed")

# Read linked files at the event boundary, while the child is still running.
live_root = root / "live-view"
p = subprocess.Popen([cli, "solve", str(root / "tiny.cas.h5"), "--backend", "cpu", "--iterations", "1", "--data", str(root / "tiny.dat.h5"),
                      "--out", str(root / "live-result"), "--progress", "jsonl", "--live-view", str(live_root), "--live-every", "1"],
                     stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
events = []
for line in p.stdout:
    event = json.loads(line)
    events.append(event)
    if event["event"] == "snapshot":
        manifest_file = pathlib.Path(event["manifest"])
        view = json.loads(manifest_file.read_text())
        assert view["v"] == 1 and len(view["zones"]) == 1 and not view["slices"]
        for zone in view["zones"]:
            data = (manifest_file.parent / zone["file"]).read_bytes()
            assert data[:4] == b"KUTV" and len(data) >= 32
        assert event["snapshot_seconds"] > 0 and event["iter"] >= 1
assert p.wait(timeout=10) == 0, p.stderr.read()
p.stdin.close()
snapshots = [e for e in events if e["event"] == "snapshot"]
assert snapshots and [s["sequence"] for s in snapshots] == list(range(1, len(snapshots) + 1))
assert events[0]["event"] == "setup" and events[-1]["event"] == "done"
config = events[0]["live_view"]
assert config["every"] == 1 and config["unit"] == "iterations" and config["budget_fraction"] == .05
pointer = json.loads((live_root / "live/manifest.json").read_text())
assert pointer["sequence"] == snapshots[-1]["sequence"] and pointer["run_id"] == config["run_id"]
assert (live_root / "live" / pointer["manifest"]).resolve() == pathlib.Path(snapshots[-1]["manifest"]).resolve()
summary = json.loads((root / "live-result/summary.json").read_text())
assert summary["live_view"] == events[-1]["live_view"] and summary["live_view"]["snapshots"] == len(snapshots)
assert math.isclose(summary["live_view"]["snapshot_time_seconds"], sum(s["snapshot_seconds"] for s in snapshots), rel_tol=1e-12)
# A finished live surface uses exactly the same fields and geometry as export.
command("export-view", root / "live-result")
final_manifest = json.loads((root / "live-result/view/manifest.json").read_text())
assert view == final_manifest and view["reference"]
for zone in view["zones"]:
    assert (manifest_file.parent / zone["file"]).read_bytes() == (root / "live-result/view" / zone["file"]).read_bytes()
# Closing stdin after observing a live update still closes ordinary run outputs.
stopping = subprocess.Popen([cli, "solve", str(root / "tiny.cas.h5"), "--backend", "cpu", "--iterations", "1000000",
                             "--out", str(root / "live-stopped"), "--progress", "jsonl", "--live-view", str(live_root), "--live-every", "1"],
                            stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
prefix = []
while not prefix or prefix[-1]["event"] != "snapshot":
    line = stopping.stdout.readline()
    assert line, stopping.stderr.read()
    prefix.append(json.loads(line))
stopping.stdin.close()
stopped_events = prefix + finish(stopping)
assert stopped_events[-1]["event"] == "done" and stopped_events[-1]["stop_reason"] == "stopped"
assert stopped_events[-1]["iterations"] < 1000000
assert stopped_events[0]["live_view"]["run_id"] != config["run_id"]
assert pathlib.Path(prefix[-1]["manifest"]).is_file()
for arguments in [("--live-every", "1"), ("--live-view", str(live_root), "--live-every", "0"),
                  ("--live-view", str(live_root), "--live-every", "1.5"),
                  ("--live-view", str(live_root), "--live-every", "99999999999999")]:
    failure = command("solve", root / "tiny.cas.h5", "--progress", "jsonl", *arguments, code=2)
    assert json.loads(failure.stdout)["event"] == "error"

# A requested cadence between monitor refreshes must never export the solver's
# pressure-correction scratch buffer. N=7 first becomes eligible at monitor 20.
refreshed = subprocess.Popen([cli, "solve", str(root / "tiny.cas.h5"), "--backend", "cpu", "--iterations", "21",
                              "--out", str(root / "live-refreshed"), "--progress", "jsonl", "--live-view", str(live_root), "--live-every", "7"],
                             stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
refreshed_events = finish(refreshed)
refreshed.stdin.close()
assert [e["iter"] for e in refreshed_events if e["event"] == "snapshot"] == [20]

terminal = subprocess.Popen([cli, "solve", str(root / "tiny.cas.h5"), "--backend", "cpu", "--iterations", "7",
                             "--out", str(root / "live-terminal"), "--progress", "jsonl", "--live-view", str(live_root), "--live-every", "7"],
                            stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
terminal_events = finish(terminal)
terminal.stdin.close()
assert [e["iter"] for e in terminal_events if e["event"] == "snapshot"] == [7]
print("Live CLI: publication, paths, run identity, sequence, timing, refreshed pressure states, stop handling and option validation passed")
