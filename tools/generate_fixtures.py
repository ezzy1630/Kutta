# /// script
# requires-python = "==3.14.*"
# dependencies = ["gmsh==4.15.2", "numpy==2.5.4", "scipy==1.18.1", "h5py==3.16.0"]
# ///
"""Generate the meshes and cases the slow tests read from the runs root.

  uv run tools/generate_fixtures.py [--runs DIR] [--check] [--force] [--record] [--only NAME]

The runs root is --runs, else $KUTTA_RUNS, else <checkout>/runs (gitignored);
tests find it the same way (tests/generated_fixtures.hpp). Generated files stay
out of git: they are large, and most derive from restricted reference data in
validation/data. validation/generated_fixtures.tsv records each file's size and
SHA-256, so a regenerated file is known to be the recorded one.

Missing fixtures are generated into a scratch directory inside the root and
moved into place only once their hash matches the ledger; a mismatch is left
in the scratch directory and reported. Existing fixtures are verified, not
rewritten (--force regenerates them). --check only verifies. --record writes
the ledger from the generated files, for a deliberate generator change.

The pinned interpreter and packages are part of the recorded identity: gmsh's
meshes depend on its version.
"""
import argparse
import hashlib
import os
import pathlib
import shutil
import subprocess
import sys
import tempfile

ROOT = pathlib.Path(__file__).resolve().parent.parent
LEDGER = ROOT / "validation/generated_fixtures.tsv"
DATA = ROOT / "validation/data"
ZERIHAN = DATA / "zerihan_tyrrell026_ground_effect/tyrrell026_single_element_coords_ref_incidence_1deg.csv"
S1223 = DATA / "uiuc_lsat_s1223/lsat_vol1_COORD01_s1223.dap.txt"
T3A = DATA / "ercoftac_t3_flat_plate/raw/t3ay.dat"


def wing(case, flag, point, name):
    # The arguments tools/fsae_aero/run.py passes, so its meshes are these.
    return ["tools/fsae_aero/wing2d.py", case, flag, point, "--level", "2", "--out", name, "--far", "1.0"]


# name (relative to the runs root) -> (generator arguments, inputs it reads).
# `{out}` is the output path; `{elbow}` the Fluent mixing elbow case.
FIXTURES = {
    "fsae_aero/mesh/ahmed35_L1.msh":
        (["tools/fsae_aero/ahmed.py", "--slant", "35", "--level", "1", "--out", "{out}"], []),
    "fsae_aero/mesh/ahmed35_L1_first15um.msh":
        (["tools/fsae_aero/ahmed.py", "--slant", "35", "--level", "1", "--out", "{out}", "--first", "1.5e-05"], []),
    "fsae_aero/mesh/s1223_a4.76_L2.msh": (wing("s1223", "--alpha", "4.76", "{out}"), [S1223]),
    "fsae_aero/mesh/zerihan_hfree_L2.msh": (wing("zerihan", "--hc", "free", "{out}"), [ZERIHAN]),
    "fsae_aero/mesh/zerihan_h0.224_L2.msh": (wing("zerihan", "--hc", "0.224", "{out}"), [ZERIHAN]),
    "fsae_aero/mesh/zerihan_h0.090_L2.msh": (wing("zerihan", "--hc", "0.090", "{out}"), [ZERIHAN]),
    # flatplate.py names the file itself: <out dir>/t3a_L1.msh.
    "transition/flatplate/t3a_L1.msh":
        (["tools/transition/flatplate.py", "mesh", "t3a", "--level", "1", "--out", "{dir}"], [T3A]),
    # The kutta CLI's 2D transition case (app/scripts/capture_fixtures.py --wing).
    "wing/naca4412_ground.cas":
        (["tools/fsae_aero/naca_ground_case.py", "--elbow", "{elbow}", "--out", "{out}"], ["{elbow}"]),
}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_ledger():
    if not LEDGER.exists():  # before the first --record
        return {name: (0, "unrecorded") for name in FIXTURES}
    rows = LEDGER.read_text().splitlines()
    if not rows or rows[0] != "path\tbytes\tsha256":
        raise SystemExit(f"{LEDGER}: invalid header")
    ledger = {}
    for row in rows[1:]:
        path, size, digest = row.split("\t")
        ledger[path] = (int(size), digest)
    if ledger.keys() != FIXTURES.keys():
        raise SystemExit(f"{LEDGER} and FIXTURES list different files: {sorted(ledger.keys() ^ FIXTURES.keys())}")
    return ledger


def write_ledger(found):
    rows = ["path\tbytes\tsha256"] + [f"{name}\t{size}\t{digest}" for name, (size, digest) in sorted(found.items())]
    LEDGER.write_text("\n".join(rows) + "\n")


def identity(path):
    return path.stat().st_size, sha256(path)


def generate(name, scratch, elbow):
    """Run the generator into `scratch`; returns the generated file."""
    command, inputs = FIXTURES[name]
    out = scratch / pathlib.Path(name).name
    missing = [p for p in (pathlib.Path(str(i).format(elbow=elbow)) for i in inputs) if not p.exists()]
    if missing:
        raise RuntimeError("missing input " + ", ".join(map(str, missing)) + " (restricted data: see its SOURCE.md)")
    args = [a.format(out=out, dir=scratch, elbow=elbow) for a in command]
    subprocess.run([sys.executable, str(ROOT / args[0])] + args[1:], check=True, cwd=ROOT, stdout=subprocess.DEVNULL)
    if not out.exists():
        raise RuntimeError(f"{args[0]} did not write {out.name}")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs", type=pathlib.Path, default=pathlib.Path(os.environ.get("KUTTA_RUNS", ROOT / "runs")))
    ap.add_argument("--fluent-examples", type=pathlib.Path,
                    default=pathlib.Path(os.environ.get("KUTTA_FLUENT_EXAMPLES", ROOT / "validation/fluent_examples")))
    ap.add_argument("--only", action="append", choices=sorted(FIXTURES), help="limit to these fixtures")
    mode = ap.add_mutually_exclusive_group()
    mode.add_argument("--check", action="store_true", help="verify existing fixtures; write nothing")
    mode.add_argument("--record", action="store_true", help="regenerate and record the ledger")
    ap.add_argument("--force", action="store_true", help="regenerate fixtures that already exist")
    a = ap.parse_args()
    ledger = read_ledger()
    names = a.only or sorted(FIXTURES)
    elbow = a.fluent_examples / "mixing_elbow/mixing_elbow.cas.h5"
    runs = a.runs.resolve()
    failed, recorded = [], dict(ledger)
    scratch = None
    try:
        for name in names:
            target = runs / name
            found = None
            if target.exists() and not (a.force or a.record):
                found = identity(target)
                status = "ok" if found == ledger[name] else "DIFFERS from the ledger (--force regenerates it)"
            elif a.check:
                status = "missing"
            else:
                if scratch is None:
                    runs.mkdir(parents=True, exist_ok=True)
                    scratch = pathlib.Path(tempfile.mkdtemp(prefix=".generate-", dir=runs))
                work = scratch / name.replace("/", "_")
                work.mkdir()
                try:
                    out = generate(name, work, elbow)
                except (RuntimeError, subprocess.CalledProcessError) as e:
                    shutil.rmtree(work)
                    failed.append(name)
                    print(f"{name}: FAILED: {e}", flush=True)
                    continue
                found = identity(out)
                if found == ledger[name] or a.record:
                    target.parent.mkdir(parents=True, exist_ok=True)
                    for f in work.iterdir():  # the file and its .json sidecar
                        if f.name.startswith(out.name):
                            os.replace(f, target.parent / f.name)
                    shutil.rmtree(work)
                    status = "generated" if found == ledger[name] else "generated, recorded"
                else:
                    # Kept for inspection; the target stays absent, so tests skip.
                    status = f"generated but DIFFERS from the ledger; kept at {out}"
            if found and a.record:
                recorded[name] = found
            if status not in ("ok", "generated", "generated, recorded"):
                failed.append(name)
            detail = f" ({found[0]} bytes, sha256 {found[1]})" if found else ""
            print(f"{name}: {status}{detail}", flush=True)
    finally:
        if scratch is not None and not any(scratch.iterdir()):
            scratch.rmdir()
    if a.record:
        write_ledger(recorded)
        print(f"recorded {LEDGER.relative_to(ROOT)}")
    print(f"{len(names) - len(failed)} of {len(names)} fixtures ok under {runs}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
