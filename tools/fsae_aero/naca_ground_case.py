#!/usr/bin/env python3
"""A 2D FSAE ground-effect wing as a Fluent text case (.cas) the kutta CLI opens.

The geometry is public: an inverted NACA 4-digit section (analytic, closed
trailing edge) above a moving ground, meshed by wing2d.build. The settings are
the Ansys mixing-elbow example's (MIT) Fluent 22.2 variables with the case
edited the way a user would in Fluent: 2D, energy off, SST with the
Intermittency Transition Model, air at the given speed, a moving ground,
drag and downforce report definitions on the wing, and the reference values
of the wing. Nothing here comes from restricted validation data.

  build-py/bin/python tools/fsae_aero/naca_ground_case.py \
      --elbow $KUTTA_FLUENT_EXAMPLES/mixing_elbow/mixing_elbow.cas.h5 --out runs/wing/naca4412_ground.cas

Needs gmsh, numpy, scipy and h5py (to read the elbow's settings).
"""
import argparse
import math
import pathlib
import re
import sys
import tempfile

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import fluent_msh  # noqa: E402
import wing2d  # noqa: E402


def naca4(code, n=200):
    """Lower TE -> LE -> upper TE of a NACA 4-digit section, unit chord, closed TE."""
    m, p, t = int(code[0]) / 100, int(code[1]) / 10, int(code[2:]) / 100
    beta = np.linspace(0, math.pi, n)
    x = 0.5 * (1 - np.cos(beta))
    yt = 5 * t * (0.2969 * np.sqrt(x) - 0.1260 * x - 0.3516 * x**2 + 0.2843 * x**3 - 0.1036 * x**4)
    front = x < p
    yc = np.where(front, m / max(p, 1e-9)**2 * (2 * p * x - x**2), m / (1 - p)**2 * ((1 - 2 * p) + 2 * p * x - x**2))
    dyc = np.where(front, 2 * m / max(p, 1e-9)**2 * (p - x), 2 * m / (1 - p)**2 * (p - x))
    th = np.arctan(dyc)
    upper = np.stack([x - yt * np.sin(th), yc + yt * np.cos(th)], axis=1)
    lower = np.stack([x + yt * np.sin(th), yc - yt * np.cos(th)], axis=1)
    return np.concatenate([lower[::-1], upper[1:]])


def balanced(text, start):
    """End index (exclusive) of the parenthesised form starting at `start`."""
    depth, i = 0, start
    while True:
        c = text[i]
        if c == '"':
            i = text.index('"', i + 1)
        elif c == "(":
            depth += 1
        elif c == ")":
            depth -= 1
            if depth == 0:
                return i + 1
        i += 1


def set_variable(rp, key, value):
    """Replace the top-level `(key value)` entry of the Rampant Variables list."""
    i = rp.find("\n(" + key + " ")
    if i < 0:
        raise SystemExit("elbow variables have no " + key)
    return rp[:i + 1] + f"({key} {value})" + rp[balanced(rp, i + 1):]


def set_config(rp, key, value):
    pattern = re.escape(f"({key} . ") + r"[^)]*\)"
    rp, n = re.subn(pattern, f"({key} . {value})", rp, count=1)
    if n != 1:
        raise SystemExit("elbow case-config has no " + key)
    return rp


def thread(threads, header):
    """The alist body of the elbow's zone section that starts with `header`."""
    i = threads.find(header)
    if i < 0:
        raise SystemExit("elbow has no zone " + header)
    body = threads.index(")(", i) + 1
    return threads[body:balanced(threads, body)]


def set_entry(alist, key, value):
    pattern = r"\n \(" + re.escape(key) + r" [^\n]*\)\n"
    # Atoms are dotted pairs, (k . 1); profiles are lists, (k (constant . 1) ...).
    entry = f"({key} {value})" if value.startswith("(") else f"({key} . {value})"
    alist, n = re.subn(pattern, lambda _: f"\n {entry}\n", alist, count=1)
    if n != 1:
        raise SystemExit("zone variables have no " + key)
    return alist


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--elbow", required=True, help="mixing_elbow.cas.h5 (settings source)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--naca", default="4412")
    ap.add_argument("--alpha", type=float, default=3.0, help="nose-down incidence, deg")
    ap.add_argument("--hc", type=float, default=0.2, help="ground clearance of the lowest point over chord")
    ap.add_argument("--chord", type=float, default=0.25, help="m")
    ap.add_argument("--speed", type=float, default=15.0, help="m/s; the ground moves at this speed")
    ap.add_argument("--intensity", type=float, default=0.01)
    ap.add_argument("--viscosity-ratio", type=float, default=10.0)
    ap.add_argument("--level", type=int, default=1, help="wing2d mesh level")
    ap.add_argument("--far", type=float, default=0.5, help="scale wing2d's far-field distances")
    a = ap.parse_args()
    import h5py

    settings = h5py.File(a.elbow, "r")["settings"]
    rp = settings["Rampant Variables"][0].decode()
    threads = settings["Thread Variables"][0].decode()
    version = settings["Version"][0].decode()

    # Inverted (camber down) and pitched nose down about the quarter chord. The
    # mirror reverses the winding; wing2d.build needs lower (suction) TE -> LE
    # -> upper TE, clockwise, so the boundary layer marches into the fluid.
    contour = (naca4(a.naca) * np.array([1.0, -1.0]))[::-1]
    contour = wing2d.rotate(contour, a.alpha, np.array([0.25, 0.0]))
    ground_y = float(contour[:, 1].min()) - a.hc
    chord = a.chord
    write = fluent_msh.write
    with tempfile.TemporaryDirectory() as tmp:
        msh = pathlib.Path(tmp) / "wing.msh"
        fluent_msh.write = lambda path, points, *rest, **kw: write(path, points * chord, *rest, **kw)
        meta = wing2d.build(contour, ground_y, a.level, msh, True,
                            {"case": f"naca{a.naca} inverted, alpha {a.alpha} deg, h/c {a.hc}"}, far=a.far)
        mesh = msh.read_text()
    mesh = "\n".join(l for l in mesh.splitlines() if not l.startswith(("(0 ", "(2 ", "(45 "))) + "\n"

    for key, value in [("rf-energy?", "#f"), ("rp-3d?", "#f"), ("rp-trans-one-eqn?", "#t"),
                       ("rp-transition-model", "gamma-transport-eqn")]:
        rp = set_config(rp, key, value)
    area = chord * 1.0
    for key, value in [("rp-trans-one-eqn?", "#t"), ("rp-transition-model", "gamma-transport-eqn"),
                       ("reference-area", f"{area:.17g}"), ("reference-depth", "1"), ("reference-length", f"{chord:.17g}"),
                       ("reference-velocity", f"{a.speed:.17g}"), ("reference-density", "1.225"),
                       # Gradients on planar 2D faces need no warped-face correction.
                       ("recon/cell-lsf?", "#f"),
                       ("monitor/report-definitions",
                        '(((name "drag") (report-definition force-drag (zone-names "wing") (force-vector 1 0 0) (scaled? #t)))'
                        ' ((name "downforce") (report-definition force-lift (zone-names "wing") (force-vector 0 -1 0) (scaled? #t))))'),
                       ("residuals/settings",
                        "(" + " ".join(f"({eq} #t 0. #t 0.0001 0. 0.05)" for eq in
                                       ("continuity", "x-velocity", "y-velocity", "k", "omega", "intermit")) + ")")]:
        rp = set_variable(rp, key, value)

    turbulence = [("ke-spec", "2"), ("turb-intensity", f"{a.intensity:.17g}"), ("turb-viscosity-ratio", f"{a.viscosity_ratio:.17g}"),
                  ("intermit", "(constant . 1) (profile \"\" \"\")")]
    inlet = thread(threads, "(39 (31 velocity-inlet cold-inlet")
    for key, value in [("velocity-spec", "2"), ("vmag", f"(constant . {a.speed:.17g}) (profile \"\" \"\")"), *turbulence]:
        inlet = set_entry(inlet, key, value)
    outlet = thread(threads, "(39 (32 pressure-outlet outlet")
    for key, value in [("p", "(constant . 0) (profile \"\" \"\")"), *turbulence]:
        outlet = set_entry(outlet, key, value)
    wall = thread(threads, "(39 (34 wall wall-elbow")
    ground = wall
    for key, value in [("motion-bc", "1"), ("relative?", "#f"), ("rotating?", "#f"),
                       ("vmag", f"{a.speed:.17g}"), ("ni", "1"), ("nj", "0"), ("nk", "0"), ("components?", "#f")]:
        ground = set_entry(ground, key, value)
    zones = [(2, "fluid", "fluid", thread(threads, "(39 (87 fluid elbow-fluid")),
             (3, "interior", "default-interior", thread(threads, "(39 (89 interior")),
             (4, "wall", "wing", wall), (5, "velocity-inlet", "inlet", inlet), (6, "pressure-outlet", "outlet", outlet),
             (7, "symmetry", "top", thread(threads, "(39 (29 symmetry")), (8, "wall", "ground", ground)]
    out = pathlib.Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as f:
        f.write(f'(0 "fluent{version}.0  generated by tools/fsae_aero/naca_ground_case.py")\n')
        f.write(f'(0 "Inverted NACA {a.naca}, {a.alpha} deg nose down, h/c {a.hc}, chord {chord} m, {a.speed} m/s, moving ground;'
                f' settings from the Ansys mixing elbow example, MIT licence")\n')
        f.write("(2 2)\n")
        f.write(rp.strip() + "\n")
        f.write(mesh)
        for zid, kind, name, alist in zones:
            f.write(f"(39 ({zid} {kind} {name} 1){alist})\n")
    print(f"wrote {out}: {meta['cells']} cells, ground y {ground_y * chord:.5f} m")


if __name__ == "__main__":
    main()
