"""Capture the app's mixing-elbow fixtures from a real `kutta` CLI.

Every file under app/fixtures/cases/mixing_elbow is real CLI output for the
Ansys mixing-elbow example (MIT), with local absolute paths replaced by the
bare file names so no machine path is committed.

Usage:
  python3 -I app/scripts/capture_elbow_fixtures.py --kutta <kutta binary> \
      --elbow-dir <dir with mixing_elbow.cas.h5 and mixing_elbow.dat.h5>
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
OUT = APP / "fixtures" / "cases" / "mixing_elbow"
SLICE = "0,0,0.025,0,0,1"  # mid-height plane through both inlets and the elbow
SNAPSHOT_ITERATIONS = (40, 120, 300)


def solve(kutta, case, data, iterations, out_dir, progress_path):
    cmd = [kutta, "solve", str(case), "--approximate", "wfgc", "--iterations", str(iterations),
           "--out", str(out_dir), "--progress", "jsonl"]
    if data:
        cmd[3:3] = ["--data", str(data)]
    with subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE) as proc:
        lines = []
        for raw in proc.stdout:
            lines.append(raw.decode())
        proc.wait()
        if proc.returncode != 0:
            sys.exit(f"solve failed ({proc.returncode}): {proc.stderr.read().decode()}")
    if progress_path:
        progress_path.write_text("".join(lines))
    return lines


def export_view(kutta, run_dir, slices):
    cmd = [kutta, "export-view", str(run_dir)]
    for s in slices:
        cmd += ["--slice", s]
    proc = subprocess.run(cmd, capture_output=True, timeout=300)
    if proc.returncode != 0:
        sys.exit(f"export-view failed: {proc.stderr.decode()}")


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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kutta", required=True)
    ap.add_argument("--elbow-dir", required=True)
    ap.add_argument("--work", default="/tmp/kutta-elbow-capture")
    args = ap.parse_args()
    elbow = Path(args.elbow_dir).resolve()
    case, data = elbow / "mixing_elbow.cas.h5", elbow / "mixing_elbow.dat.h5"
    work = Path(args.work)
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir(parents=True)
    shutil.rmtree(OUT, ignore_errors=True)
    OUT.mkdir(parents=True)

    for name, extra in (("inspect-data.json", ["--data", str(data)]), ("inspect.json", [])):
        proc = subprocess.run([args.kutta, "inspect", str(case), *extra, "--json"],
                              capture_output=True, timeout=120)
        doc = json.loads(proc.stdout)
        (OUT / name).write_text(json.dumps(doc, separators=(",", ":")))

    # Unsupported without the WFGC opt-in: keep the CLI's own rejection.
    proc = subprocess.run([args.kutta, "solve", str(case), "--out", str(work / "rejected"),
                           "--progress", "jsonl"], capture_output=True, stdin=subprocess.PIPE, timeout=120)
    (OUT / "rejected.jsonl").write_text(proc.stdout.decode())

    variants = {"seeded": data, "cold": None}
    for variant, seed in variants.items():
        run_dir = work / f"run-{variant}"
        solve(args.kutta, case, seed, 600, run_dir, OUT / f"progress-{variant}.jsonl")
        export_view(args.kutta, run_dir, [SLICE])
        dest = OUT / f"run-{variant}"
        dest.mkdir()
        for f in ("summary.json", "forces.csv", "fluxes.csv", "residuals.csv"):
            shutil.copy(run_dir / f, dest / f)
        shutil.copytree(run_dir / "view", dest / "view")

    for it in SNAPSHOT_ITERATIONS:
        run_dir = work / f"snap-{it}"
        solve(args.kutta, case, None, it, run_dir, None)
        export_view(args.kutta, run_dir, [])
        shutil.copytree(run_dir / "view", OUT / "live" / f"iter-{it:04d}")

    devices = subprocess.run([args.kutta, "devices", "--json"], capture_output=True, timeout=60)
    (APP / "fixtures" / "devices-macos.json").write_text(devices.stdout.decode())

    replacements = [(str(case), "mixing_elbow.cas.h5"), (str(data), "mixing_elbow.dat.h5"),
                    (str(work) + "/", ""), (str(work), ".")]
    for path in OUT.rglob("*"):
        if path.suffix in (".json", ".jsonl", ".csv"):
            scrub(path, replacements)
    qcompress(OUT)
    print("captured", sum(p.stat().st_size for p in OUT.rglob("*") if p.is_file()), "bytes")


if __name__ == "__main__":
    main()
