"""Record the fixture CLI's cases from a real `kutta` CLI.

Every file under app/fixtures/cases/<case> is real CLI output, with local
absolute paths replaced by bare file names so no machine path is committed:

  mixing_elbow       Ansys mixing elbow (MIT): SST, Coupled, 3D; cold and
                     Fluent-seeded runs, live snapshots, the rejected run
  elbow_param        Ansys elbow_param: standard k-epsilon, SIMPLE, 3D
  cylinder_tutorial  Ansys adjoint cylinder tutorial: laminar, 2D
  naca4412_ground    tools/fsae_aero/naca_ground_case.py: 2D ground-effect
                     wing, SST + intermittency transition, moving ground

Usage:
  python3 -I app/scripts/capture_fixtures.py --kutta <kutta binary> \
      --examples <validation/fluent_examples> [--wing <naca4412_ground.cas>] [--cases a,b]
"""

import argparse
import json
import shutil
import struct
import subprocess
import zlib
import sys
from pathlib import Path

APP = Path(__file__).resolve().parent.parent
CASES = APP / "fixtures" / "cases"
SLICE = "0,0,0.025,0,0,1"  # elbow: mid-height plane through both inlets and the elbow

# Runs are recorded with convergence checking off; the fixture evaluates the
# app's criteria on the recorded residuals. Iterations cover each case's
# convergence at the default criteria.
SPECS = {
    "mixing_elbow": {"case": "mixing_elbow/mixing_elbow.cas.h5", "data": "mixing_elbow/mixing_elbow.dat.h5",
                     "approximate": "wfgc", "iterations": 600, "slices": [SLICE], "snapshots": (40, 120, 300)},
    "elbow_param": {"case": "mixing_elbow/elbow_param.cas.h5", "approximate": "compressibility,ke-options",
                    "iterations": 300, "slices": "default"},
    "cylinder_tutorial": {"case": "adjoint_cylinder_tutorial/cylinder_tutorial.cas", "approximate": "reverse-flow",
                          "iterations": 300},
    "naca4412_ground": {"case": None, "iterations": 600},
}


def solve(kutta, case, data, approximate, iterations, out_dir, progress_path):
    cmd = [kutta, "solve", str(case), "--iterations", str(iterations), "--out", str(out_dir), "--progress", "jsonl"]
    if approximate:
        cmd += ["--approximate", approximate]
    if data:
        cmd += ["--data", str(data)]
    with subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE) as proc:
        lines = [raw.decode() for raw in proc.stdout]
        proc.wait()
        if proc.returncode != 0:
            sys.exit(f"solve failed ({proc.returncode}): {proc.stderr.read().decode()}")
    if progress_path:
        progress_path.write_text("".join(lines))


def export_view(kutta, run_dir, slices):
    cmd = [kutta, "export-view", str(run_dir)]
    for s in slices:
        cmd += ["--slice", s]
    proc = subprocess.run(cmd, capture_output=True, timeout=300)
    if proc.returncode != 0:
        sys.exit(f"export-view failed: {proc.stderr.decode()}")


def default_slice(inspect):
    """The plane the app slices a 3D internal flow with (Session::defaultSlices):
    through the domain's centre, normal to its thinnest direction."""
    lo, hi = inspect["mesh"]["bounds"]
    ext = [b - a for a, b in zip(lo, hi)]
    thin = ext.index(min(ext))
    normal = [1 if a == thin else 0 for a in range(3)]
    return ",".join(f"{v:.9g}" for v in [(a + b) / 2 for a, b in zip(lo, hi)] + normal)


def drop_back_cap(view):
    """A 2D case's view has both extrusion caps, the same cells twice; the app
    shows only the cap at the top of the depth (SceneData makePlanar), so the
    recording keeps that one alone."""
    manifest = view / "manifest.json"
    m = json.loads(manifest.read_text())
    zmin = m["bounds"][0][2]
    keep = []
    for z in m["zones"]:
        raw = (view / z["file"]).read_bytes()
        first_z = struct.unpack_from("<f", raw, m["format"]["header_bytes"] + 8)[0]
        if z["id"] < 0 and z["type"] == "symmetry" and first_z == zmin:
            for f in (z["file"], z.get("edges", {}).get("file")):
                if f:
                    (view / f).unlink()
        else:
            keep.append(z)
    m["zones"] = keep
    manifest.write_text(json.dumps(m, separators=(",", ":")))


def scrub(path, replacements):
    text = path.read_text()
    for old, new in replacements:
        text = text.replace(old, new)
    path.write_text(text)


def qcompress(root):
    """Store bulky files in Qt's qCompress layout: uint32 BE length + zlib."""
    for path in list(root.rglob("*")):
        if path.is_file() and (path.suffix == ".bin" or path.name.startswith("inspect")):
            raw = path.read_bytes()
            path.with_name(path.name + ".qz").write_bytes(struct.pack(">I", len(raw)) + zlib.compress(raw, 9))
            path.unlink()


def capture(kutta, name, spec, case, data, work):
    out = CASES / name
    shutil.rmtree(out, ignore_errors=True)
    out.mkdir(parents=True)
    inspects = [("inspect.json", [])] + ([("inspect-data.json", ["--data", str(data)])] if data else [])
    for file, extra in inspects:
        proc = subprocess.run([kutta, "inspect", str(case), *extra, "--json"], capture_output=True, timeout=120)
        (out / file).write_text(json.dumps(json.loads(proc.stdout), separators=(",", ":")))
    inspect = json.loads((out / "inspect.json").read_text())
    slices = [default_slice(inspect)] if spec.get("slices") == "default" else spec.get("slices", [])
    planar = any(e["path"] == "general/dimension" and e["value"] == 2
                 for e in inspect["settings"])

    if spec.get("approximate"):
        # Without the opt-in the CLI rejects the case: keep its own rejection.
        proc = subprocess.run([kutta, "solve", str(case), "--out", str(work / "rejected"), "--progress", "jsonl"],
                              capture_output=True, stdin=subprocess.PIPE, timeout=120)
        (out / "rejected.jsonl").write_text(proc.stdout.decode())

    for variant, seed in (("cold", None), ("seeded", data)):
        if variant == "seeded" and not data:
            continue
        run_dir = work / f"run-{variant}"
        solve(kutta, case, seed, spec.get("approximate"), spec["iterations"], run_dir, out / f"progress-{variant}.jsonl")
        export_view(kutta, run_dir, slices)
        dest = out / f"run-{variant}"
        dest.mkdir()
        for f in ("summary.json", "forces.csv", "fluxes.csv", "residuals.csv", "loads.json"):
            shutil.copy(run_dir / f, dest / f)
        shutil.copytree(run_dir / "view", dest / "view")
        if planar:
            drop_back_cap(dest / "view")

    for it in spec.get("snapshots", ()):
        run_dir = work / f"snap-{it}"
        solve(kutta, case, None, spec.get("approximate"), it, run_dir, None)
        export_view(kutta, run_dir, [])
        shutil.copytree(run_dir / "view", out / "live" / f"iter-{it:04d}")

    replacements = [(str(case), case.name)] + ([(str(data), data.name)] if data else []) + [(str(work) + "/", ""), (str(work), ".")]
    for path in out.rglob("*"):
        if path.suffix in (".json", ".jsonl", ".csv"):
            scrub(path, replacements)
    qcompress(out)
    print(name, sum(p.stat().st_size for p in out.rglob("*") if p.is_file()), "bytes")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kutta", required=True)
    ap.add_argument("--examples", required=True, help="validation/fluent_examples")
    ap.add_argument("--wing", help="naca4412_ground.cas from tools/fsae_aero/naca_ground_case.py")
    ap.add_argument("--cases", default=",".join(SPECS))
    ap.add_argument("--work", default="/tmp/kutta-fixture-capture")
    args = ap.parse_args()
    examples = Path(args.examples).resolve()
    for name in args.cases.split(","):
        spec = SPECS[name]
        if spec["case"]:
            case = examples / spec["case"]
        elif args.wing:
            case = Path(args.wing).resolve()
        else:
            sys.exit(f"{name} needs --wing")
        data = examples / spec["data"] if spec.get("data") else None
        work = Path(args.work) / name
        shutil.rmtree(work, ignore_errors=True)
        work.mkdir(parents=True)
        capture(args.kutta, name, spec, case, data, work)
    devices = subprocess.run([args.kutta, "devices", "--json"], capture_output=True, timeout=60)
    (APP / "fixtures" / "devices-macos.json").write_text(devices.stdout.decode())


if __name__ == "__main__":
    main()
