#!/usr/bin/env python3
"""Run FSAE aero validation sweeps through build/bench_fsae_aero.

Each point gets a run directory with the solver log, history.csv, walls.csv and
result.json (RESULT line + mesh metadata + peak memory from /usr/bin/time -l).
Meshes are generated on demand into runs/fsae_aero/mesh (gitignored).

  build-py/bin/python tools/fsae_aero/run.py zerihan --points free,0.224 --level 2
  build-py/bin/python tools/fsae_aero/run.py s1223 --re 2e5 --points 0,4.76 --level 2
  build-py/bin/python tools/fsae_aero/run.py ahmed --points 35 --level 1
"""
import argparse
import json
import pathlib
import re
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
PY = sys.executable
MESHES = ROOT / "runs/fsae_aero/mesh"

# Ahmed body: U = 1, lengths in metres, Re_L = 2.78e6 (Lienhart & Becker, 40 m/s),
# full model, frontal area 0.112 m^2; inlet Tu 1% and nu_t/nu 1 as CFD11 (CFX-5 SST).
AHMED_LENGTH, AHMED_RE, AHMED_AREA = 1.044, 2.78e6, 0.389 * 0.288


def mesh_for(case, point, level, far=1.0, first=None):
    if case == "ahmed":
        # point = slant angle; a non-default first prism height is in the name
        suffix = f"_first{first * 1e6:g}um" if first else ""
        path = MESHES / f"ahmed{point}_L{level:g}{suffix}.msh"
        cmd = [PY, str(ROOT / "tools/fsae_aero/ahmed.py"), "--slant", point, "--level", f"{level:g}", "--out", str(path)]
        cmd += ["--first", str(first)] if first else []
    else:
        tag = f"h{point}" if case == "zerihan" else f"a{point}"
        path = MESHES / f"{case}_{tag}_L{level:g}{f'_far{far:g}' if far != 1 else ''}.msh"
        flag = "--hc" if case == "zerihan" else "--alpha"
        cmd = [PY, str(ROOT / "tools/fsae_aero/wing2d.py"), case, flag, point, "--level", f"{level:g}", "--out", str(path),
               "--far", str(far)]
    if not path.exists():
        MESHES.mkdir(parents=True, exist_ok=True)
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
    return path


def solve(case, point, args):
    mesh = mesh_for(case, point, args.level, args.far, args.first)
    first = f"_first{args.first * 1e6:g}um" if case == "ahmed" and args.first else ""
    tag = f"{case}{'_re' + args.re if case == 's1223' else ''}_{point}_L{args.level:g}{first}_{args.backend}{args.tag}"
    out = pathlib.Path(args.out) / tag
    out.mkdir(parents=True, exist_ok=True)
    cmd = [str(ROOT / "build/bench_fsae_aero"), "--mesh", str(mesh), "--backend", args.backend, "--out", str(out),
           "--iterations", str(args.iterations), "--max-seconds", str(args.max_seconds), "--check", str(args.check),
           "--tu", str(args.tu if case != "ahmed" else 0.01), "--ratio", str(args.ratio)]
    if case == "zerihan":
        cmd += ["--re", "4.62e5", "--downforce"] + (["--moving-ground"] if point != "free" else [])
    elif case == "s1223":
        cmd += ["--re", args.re]
    else:
        cmd += ["--re", str(AHMED_RE), "--length", str(AHMED_LENGTH), "--area", str(AHMED_AREA),
                "--walls", "body", "--lift-z", "--window", "1000"]
    cmd += args.extra.split() if args.extra else []
    with open(out / "solve.log", "w") as log:
        proc = subprocess.run(["/usr/bin/time", "-l"] + cmd, stdout=log, stderr=subprocess.STDOUT)
    text = (out / "solve.log").read_text()
    m = re.search(r"^RESULT (\{.*\})$", text, re.M)
    result = json.loads(m.group(1)) if m else {"error": text[-2000:]}
    for key, pattern in (("peak_footprint_bytes", r"(\d+)\s+peak memory footprint"),
                         ("max_rss_bytes", r"(\d+)\s+maximum resident set size"),
                         ("wall_seconds", r"([\d.]+) real")):
        hit = re.search(pattern, text)
        if hit:
            result[key] = float(hit.group(1))
    result.update(case=case, point=point, level=args.level, backend=args.backend, exit=proc.returncode, tag=args.tag,
                  command=[str(pathlib.Path(c).relative_to(ROOT)) if c.startswith(str(ROOT)) else c for c in cmd],
                  mesh=json.loads(pathlib.Path(str(mesh) + ".json").read_text()))
    if case == "s1223":
        result["reynolds"] = float(args.re)
    (out / "result.json").write_text(json.dumps(result, indent=1))
    with open(pathlib.Path(args.out) / "results.jsonl", "a") as all_results:
        all_results.write(json.dumps(result) + "\n")
    r = result.get("result", {})
    print(f"{tag}: cells={r.get('cells')} it={r.get('iterations')} lift={r.get('lift')} drag={r.get('drag')} "
          f"res={r.get('max_residual')} settled={r.get('settled')} t={r.get('solve_s')} "
          f"mem={result.get('peak_footprint_bytes', 0) / 1e9:.2f}GB failure={r.get('failure', result.get('error', ''))[:200]}",
          flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case", choices=["zerihan", "s1223", "ahmed"])
    ap.add_argument("--points", required=True, help="comma list: h/c (or 'free'), alpha deg, or slant deg")
    ap.add_argument("--re", default="2e5", help="S1223 chord Reynolds number")
    ap.add_argument("--level", type=float, default=2, help="mesh level (Ahmed accepts fractions, e.g. 3.5)")
    ap.add_argument("--backend", default="metal")
    ap.add_argument("--iterations", type=int, default=6000)
    ap.add_argument("--max-seconds", type=float, default=3600)
    ap.add_argument("--check", type=int, default=50)
    ap.add_argument("--tu", type=float, default=0.001)
    ap.add_argument("--ratio", type=float, default=1.0)
    ap.add_argument("--extra", default="", help="extra bench_fsae_aero options")
    ap.add_argument("--far", type=float, default=1.0, help="2D far-field distance scale (wing2d.py --far)")
    ap.add_argument("--first", type=float, default=None, help="Ahmed first prism height, m (ahmed.py --first)")
    ap.add_argument("--tag", default="", help="suffix for the run directory")
    ap.add_argument("--out", default=str(ROOT / "runs/fsae_aero/results"))
    a = ap.parse_args()
    for point in a.points.split(","):
        solve(a.case, point, a)


if __name__ == "__main__":
    main()
