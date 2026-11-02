#!/usr/bin/env python3
"""2D wing-section meshes for the FSAE aero validation (Zerihan ground effect, S1223).

Writes an ASCII Fluent mesh (unit chord, quarter-chord pivot) with boundary zones
`wing`, `ground` (wall, or symmetry in free air), `inlet`, `outlet`, `top`.
Geometry comes from the restricted reference files under validation/data; nothing
derived from them is committed. Run with the venv that has gmsh/numpy/scipy:

  build-py/bin/python tools/fsae_aero/wing2d.py zerihan --hc 0.224 --level 2 --out runs/fsae_aero/z.msh
  build-py/bin/python tools/fsae_aero/wing2d.py s1223 --alpha 4 --level 2 --out runs/fsae_aero/s.msh

Fractional levels permit refinement within a memory limit. --ground-layers
adds a wall-resolving strip; --domain X_IN X_OUT Y_TOP changes the far boundaries
independently of local spacing. Both are opt-in and preserve the default mesh.
"""
import argparse
import csv
import json
import math
import pathlib
import sys

import numpy as np
from scipy.interpolate import CubicSpline

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import fluent_msh  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parents[2]
ZERIHAN = ROOT / "validation/data/zerihan_tyrrell026_ground_effect/tyrrell026_single_element_coords_ref_incidence_1deg.csv"
S1223 = ROOT / "validation/data/uiuc_lsat_s1223/lsat_vol1_COORD01_s1223.dap.txt"

# Domain in chords (scaled by --far). Uniform inlet and slip top/bottom:
# doubling every distance changed the freestream Cl by +0.15% (REPORT.md).
X_IN, X_OUT, Y_FAR = -30.0, 40.0, 30.0
# Recombining gmsh's triangles produced occasional inverted quads (checked by
# fluent_msh.write), so the outer region stays triangular.
MESH_OPTIONS = {"Mesh.MeshSizeExtendFromBoundary": 0, "Mesh.MeshSizeFromPoints": 0,
                "Mesh.MeshSizeFromCurvature": 0, "Mesh.Algorithm": 6, "Mesh.RecombineAll": 0,
                "Mesh.Smoothing": 5}


def rotate(p, deg, pivot):
    a = math.radians(deg)
    r = np.array([[math.cos(a), -math.sin(a)], [math.sin(a), math.cos(a)]])
    return (p - pivot) @ r.T + pivot


def zerihan_contour():
    """Lower (suction) TE -> LE -> upper TE at alpha_ref = 1 deg (true 3.45 deg)."""
    rows = list(csv.DictReader(open(ZERIHAN)))
    side = {s: np.array([[float(r["x_over_c"]), float(r["y_over_c"])] for r in rows if r["surface"] == s])
            for s in ("suction", "pressure")}
    contour = np.concatenate([side["suction"][::-1], side["pressure"][1:]])
    # The table's x is the horizontal projection; scale to unit true chord
    # (LE to trailing-edge midpoint), SOURCE.md.
    te_mid = 0.5 * (contour[0] + contour[-1])
    return contour / np.hypot(*te_mid)


def selig_contour(path=S1223):
    rows = [l.split() for l in open(path).read().splitlines()[1:] if l.strip()]
    pts = np.array([[float(a), float(b)] for a, b in rows])
    return pts[::-1]  # Selig order is upper TE -> LE -> lower TE


def vinokur(n, ds0, ds1, length):
    """n+1 points on [0, length] with end spacings ~ds0, ds1 (Vinokur 1983, two-sided)."""
    s0, s1 = length / (n * ds0), length / (n * ds1)
    b = math.sqrt(s0 * s1)
    a = math.sqrt(s0 / s1)
    xi = np.arange(n + 1) / n
    if b > 1.001:
        dy = math.asinh(b) if b < 1e6 else math.log(2 * b)
        # Solve sinh(dy)/dy = b
        for _ in range(100):
            dy -= (math.sinh(dy) / dy - b) / ((math.cosh(dy) * dy - math.sinh(dy)) / dy**2)
        u = 0.5 * (1 + np.tanh(dy * (xi - 0.5)) / math.tanh(dy / 2))
    elif b < 0.999:
        dx = math.asin(b) if b <= 1 else 1.0
        for _ in range(100):
            dx -= (math.sin(dx) / dx - b) / ((math.cos(dx) * dx - math.sin(dx)) / dx**2)
        u = 0.5 * (1 + np.tan(dx * (xi - 0.5)) / math.tan(dx / 2))
    else:
        u = xi
    t = u / (a + (1 - a) * u)
    return t * length


def resample(contour, level, closed_te_tol=1e-3):
    """Spline the contour and redistribute it; returns lower+upper points and TE info."""
    gap = np.linalg.norm(contour[0] - contour[-1])
    contour = contour.copy()
    blunt = gap > closed_te_tol
    if 0 < gap <= closed_te_tol:
        # Close a sub-0.1%c gap by shearing the last 10% of chord; the change is
        # below the coordinate tolerance of the source tables.
        mid = 0.5 * (contour[0] + contour[-1])
        for end in (0, -1):
            delta = mid - contour[end]
            w = np.clip((contour[:, 0] - 0.9) / 0.1, 0, 1)
            side = np.arange(len(contour)) < np.argmin(contour[:, 0]) if end == 0 else \
                np.arange(len(contour)) > np.argmin(contour[:, 0])
            contour[side] += (w[side] ** 2)[:, None] * delta
    s = np.concatenate([[0], np.cumsum(np.linalg.norm(np.diff(contour, axis=0), axis=1))])
    le = int(np.argmin(contour[:, 0]))
    spline = CubicSpline(s, contour, bc_type="not-a-knot")
    r = 2.0 ** (level - 2)
    n_side = int(round(120 * r))
    ds_le, ds_te = 0.0016 / r, 0.0015 / r
    lower = spline(vinokur(n_side, ds_te, ds_le, s[le]))
    upper = spline(s[le] + vinokur(n_side, ds_le, ds_te, s[-1] - s[le]))
    pts = np.concatenate([lower, upper[1:]])
    if not blunt:
        pts[-1] = pts[0]
    return pts, blunt


def preview(path, points, blocks, windows):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    fig, axes = plt.subplots(1, len(windows), figsize=(6 * len(windows), 5))
    for ax, (x0, x1, y0, y1) in zip(np.atleast_1d(axes), windows):
        for _, b in blocks:
            p = points[b][:, :, :2]
            c = p.mean(axis=1)
            keep = (c[:, 0] > x0) & (c[:, 0] < x1) & (c[:, 1] > y0) & (c[:, 1] < y1)
            p = p[keep]
            segs = np.concatenate([np.stack([p[:, k], p[:, (k + 1) % p.shape[1]]], axis=1)
                                   for k in range(p.shape[1])])
            ax.add_collection(LineCollection(segs, linewidths=0.2, colors="k"))
        ax.set_xlim(x0, x1)
        ax.set_ylim(y0, y1)
        ax.set_aspect("equal")
    fig.savefig(path, dpi=200)
    plt.close(fig)


def oblock(loop, ds, ratio, thickness, top_size):
    """Structured near-wall O-block marched outward from a closed wing loop.

    `loop` runs clockwise (lower TE -> LE -> upper TE -> base), so the left
    normal points into the fluid. Each layer re-evaluates normals from the
    previous one; convex corners (the blunt base) fan out by construction.
    Growth is geometric from `ds` and capped at `top_size`, so the last layer
    is nearly isotropic where the unstructured mesh takes over.
    """
    layers = [loop]
    h, total = ds, 0.0
    while total + h <= thickness:
        p = layers[-1]
        # Bisector of the two adjacent edge normals, independent of the edge
        # lengths, so unequal spacing at a corner cannot tilt the line.
        e = np.roll(p, -1, axis=0) - p
        e /= np.linalg.norm(e, axis=1)[:, None]
        en = np.stack([-e[:, 1], e[:, 0]], axis=1)
        normal = en + np.roll(en, 1, axis=0)
        normal /= np.linalg.norm(normal, axis=1)[:, None]
        q = p + h * normal
        total += h
        # Tangential smoothing grows with distance from the wall: orthogonal
        # first cells, then a smooth outer front around the corners.
        w = 0.25 * min(1.0, total / thickness)
        for _ in range(4):
            mid = 0.5 * (np.roll(q, 1, axis=0) + np.roll(q, -1, axis=0))
            q = q + w * (mid - q)
        layers.append(q)
        h = min(h * ratio, top_size)
    return np.array(layers)


def signed_area(points, cells):
    x, y = points[cells, 0], points[cells, 1]
    return 0.5 * np.sum(x * np.roll(y, -1, axis=1) - np.roll(x, -1, axis=1) * y, axis=1)


def build(contour, ground_y, level, out, ground_wall, meta, png=None, far=1.0, ground_layers=False, domain=None,
          wake_refinement=1.0, trailing_edge_refinement=1, wall_refinement=1.0):
    import gmsh
    pts, blunt = resample(contour, level + math.log2(wall_refinement))
    r = 2.0 ** (level - 2)
    ds_wall = 2e-5 / (r * wall_refinement)
    ratio = 1.2 ** (1.0 / (r * wall_refinement))
    surf = 0.004 / r
    n = len(pts) - (0 if blunt else 1)
    # Clockwise loop: lower TE -> LE -> upper TE, then down the blunt base.
    loop = pts[:n]
    if blunt:
        # Uniform base spacing close to the trailing-edge surface spacing.
        nb = max(4, int(round(np.linalg.norm(loop[0] - loop[-1]) / (0.0012 / (r * wall_refinement)))))
        base = loop[-1] + (loop[0] - loop[-1]) * (np.arange(1, nb) / nb)[:, None]
        loop = np.concatenate([loop, base])
    if trailing_edge_refinement != 1:
        # Split existing straight wall faces, preserving the polygon exactly.
        # This isolates the last 15% chord/base resolution from geometry changes.
        te_mid = 0.5 * (pts[0] + pts[-1])
        leading_edge = pts[np.argmax(np.linalg.norm(pts - te_mid, axis=1))]
        chord = te_mid - leading_edge
        chord_length_squared = chord @ chord
        refined = []
        for p, q in zip(loop, np.roll(loop, -1, axis=0)):
            x0 = (p - leading_edge) @ chord / chord_length_squared
            x1 = (q - leading_edge) @ chord / chord_length_squared
            count = trailing_edge_refinement if min(x0, x1) > 0.85 else 1
            refined.extend(p + (q - p) * j / count for j in range(count))
        loop = np.array(refined)
    gap = float(pts[:, 1].min() - ground_y) if ground_wall else None
    thickness = 0.02 if gap is None else min(0.02, 0.3 * gap)
    block = oblock(loop, ds_wall, ratio, thickness, surf)
    nl, nloop = block.shape[0], block.shape[1]
    hole = block[-1]

    gmsh.initialize(["-nt", "1"])
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.option.setNumber("General.NumThreads", 1)
    gmsh.model.add("wing")
    geo = gmsh.model.geo
    tags = [geo.addPoint(x, y, 0) for x, y in hole]
    ring = [geo.addLine(tags[a], tags[(a + 1) % nloop]) for a in range(nloop)]
    x_in, x_out, y_far = domain if domain else (X_IN * far, X_OUT * far, Y_FAR * far)
    corner_out, corner_in = geo.addPoint(x_out, y_far, 0), geo.addPoint(x_in, y_far, 0)
    # Optional ground strip: share every node of its top with gmsh,
    # avoiding the terminating layers/slivers of a local ground block.
    ground_offsets = [0.0]
    if ground_layers and ground_wall:
        dy = ds_wall
        while ground_offsets[-1] + dy <= thickness:
            ground_offsets.append(ground_offsets[-1] + dy)
            dy = min(dy * ratio, surf)
    def ground_blend(x):
        t = np.clip(np.maximum(-1.0 - x, x - 4.0) / 3.0, 0, 1)
        return t * t * (3 - 2 * t)

    if len(ground_offsets) > 1:
        # Grow the whole strip in the far field, rather than carrying very
        # thin cells under metre-long triangles. This is flat through the
        # wing and near wake, where the ground boundary layer matters.
        xs = [x_in]
        while xs[-1] < x_out:
            distance = max(-0.5 - xs[-1], xs[-1] - 3.0, 0)
            dx = min(2.0, min(0.01, gap / 12) / r + 0.08 * distance)
            xs.append(min(x_out, xs[-1] + dx))
        gx = np.array(xs)
        height = ground_offsets[-1] + (1.0 - ground_offsets[-1]) * ground_blend(gx)
        ground_tags = [geo.addPoint(x, ground_y + h, 0) for x, h in zip(gx, height)]
        ground = [geo.addLine(a, b) for a, b in zip(ground_tags[:-1], ground_tags[1:])]
        start, end = ground_tags[0], ground_tags[-1]
    else:
        start, end = geo.addPoint(x_in, ground_y, 0), geo.addPoint(x_out, ground_y, 0)
        ground = [geo.addLine(start, end)]
    outlet = geo.addLine(end, corner_out)
    top = geo.addLine(corner_out, corner_in)
    inlet = geo.addLine(corner_in, start)
    outer = geo.addCurveLoop(ground + [outlet, top, inlet])
    inner = geo.addCurveLoop(ring)
    fluid = geo.addPlaneSurface([outer, inner])
    geo.synchronize()
    m = gmsh.model.mesh
    for c in ring:
        m.setTransfiniteCurve(c, 2)
    if len(ground_offsets) > 1:
        for c in ground:
            m.setTransfiniteCurve(c, 2)
    field = m.field
    # Background sizes: near-body, wake and (in ground effect) the gap.
    dist = field.add("Distance")
    field.setNumbers(dist, "CurvesList", ring)
    field.setNumber(dist, "Sampling", 4)
    thr = field.add("Threshold")
    field.setNumber(thr, "InField", dist)
    field.setNumber(thr, "SizeMin", surf)
    field.setNumber(thr, "SizeMax", 2.0)
    field.setNumber(thr, "DistMin", 0.0)
    field.setNumber(thr, "DistMax", 25.0)
    field.setNumber(thr, "StopAtDistMax", 0)
    sizes = [thr]
    te = 0.5 * (pts[0] + pts[n - 1])
    for x1, half, size in ((3.0, 0.15, 0.008), (10.0, 0.6, 0.03)):
        box = field.add("Box")
        field.setNumber(box, "VIn", size / r)
        field.setNumber(box, "VOut", 2.0)
        field.setNumber(box, "XMin", float(pts[:, 0].min()) - 0.05)
        field.setNumber(box, "XMax", te[0] + x1)
        field.setNumber(box, "YMin", te[1] - half - (0.1 * x1 if ground_wall else 0))
        field.setNumber(box, "YMax", te[1] + half + 0.1 * x1)
        field.setNumber(box, "Thickness", 12.0)  # ~15% growth outside the box
        sizes.append(box)
    if wake_refinement != 1:
        box = field.add("Box")
        field.setNumber(box, "VIn", 0.008 / (r * wake_refinement))
        field.setNumber(box, "VOut", 2.0)
        field.setNumber(box, "XMin", float(te[0]))
        field.setNumber(box, "XMax", float(te[0]) + 3)
        # Cover the deflected free-air wake on either side of the trailing edge.
        field.setNumber(box, "YMin", float(te[1]) - 0.45)
        field.setNumber(box, "YMax", float(te[1]) + 0.45)
        field.setNumber(box, "Thickness", 12.0)
        sizes.append(box)
    if ground_wall:
        box = field.add("Box")
        field.setNumber(box, "VIn", min(0.01, gap / 12) / r)
        field.setNumber(box, "VOut", 2.0)
        field.setNumber(box, "XMin", -0.5)
        field.setNumber(box, "XMax", 3.0)
        field.setNumber(box, "YMin", ground_y)
        field.setNumber(box, "YMax", float(pts[:, 1].max()))
        field.setNumber(box, "Thickness", 12.0)
        sizes.append(box)
    fmin = field.add("Min")
    field.setNumbers(fmin, "FieldsList", sizes)
    field.setAsBackgroundMesh(fmin)
    for key, value in MESH_OPTIONS.items():
        gmsh.option.setNumber(key, value)
    groups = {"inlet": [inlet], "outlet": [outlet], "top": [top], "ground": ground}
    for name, curves in groups.items():
        gmsh.model.addPhysicalGroup(1, curves, name=name)
    gmsh.model.addPhysicalGroup(2, [fluid], name="fluid")
    m.generate(2)
    node_tags, coords, _ = m.getNodes()
    gpoints = np.zeros((int(node_tags.max()) + 1, 3))
    gpoints[node_tags.astype(np.int64)] = coords.reshape(-1, 3)
    types, _, enodes = m.getElements(2)
    gblocks = [(int(t), nodes.astype(np.int64).reshape(-1, {2: 3, 3: 4}[int(t)])) for t, nodes in zip(types, enodes)]
    gfaces = {}
    for name, curves in groups.items():
        gfaces[name] = np.concatenate([m.getElements(1, c)[2][0].astype(np.int64).reshape(-1, 2) for c in curves])
    # gmsh nodes on the structured interfaces are the blocks' last layers, bit for bit.
    node_of = lambda t: int(m.getNodes(0, t)[0][0])
    ring_nodes = [node_of(t) for t in tags]
    gmsh.finalize()
    if not np.array_equal(gpoints[ring_nodes, :2], hole):
        raise SystemExit("gmsh moved structured-block interface nodes")

    points = [gpoints]
    count = len(gpoints)

    def attach(grid, shared):
        # Nodes shared with gmsh keep gmsh's ids (shared[k, j] >= 0); the rest are new.
        nonlocal count
        ids = shared.copy()
        new = ids < 0
        ids[new] = count + np.arange(np.count_nonzero(new))
        count += np.count_nonzero(new)
        points.append(np.concatenate([grid[new], np.zeros((np.count_nonzero(new), 1))], axis=1))
        return ids

    gsign = np.sign(np.sum(np.sign(signed_area(gpoints, gblocks[0][1]))))

    def quads_of(ids, closed):
        a = np.arange(ids.shape[1] if closed else ids.shape[1] - 1)
        b = (a + 1) % ids.shape[1]
        q = np.stack([ids[:-1, a], ids[:-1, b], ids[1:, b], ids[1:, a]], axis=2).reshape(-1, 4)
        # Match gmsh's cell orientation so the converter sees one consistent sign.
        return q if np.sign(np.sum(np.sign(signed_area(np.concatenate(points), q)))) == gsign else q[:, ::-1]

    shared = -np.ones(block.shape[:2], dtype=np.int64)
    shared[-1] = ring_nodes
    ids = attach(block, shared)
    blocks = gblocks + [(3, quads_of(ids, True))]
    a = np.arange(nloop)
    boundaries = [("wing", "wall", np.stack([ids[0, a], ids[0, (a + 1) % nloop]], axis=1))]
    if len(ground_offsets) > 1:
        interface = np.unique(gfaces["ground"])
        interface = interface[np.argsort(gpoints[interface, 0])]
        strip = np.empty((len(ground_offsets), len(interface), 2))
        strip[:, :, 0] = gpoints[interface, 0]
        # Retain the wall resolution through the wing and near wake. Farther
        # out, smoothly relax to uniform layers in a thicker strip so the
        # far-field cells do not inherit a 1e-5c height.
        x = gpoints[interface, 0]
        blend = ground_blend(x)
        height = gpoints[interface, 1] - ground_y
        uniform = np.linspace(0, 1, len(ground_offsets))
        fractions = np.array(ground_offsets) / ground_offsets[-1]
        offsets = (fractions[:, None] * (1 - blend) + uniform[:, None] * blend) * height
        strip[:, :, 1] = ground_y + offsets
        shared = -np.ones(strip.shape[:2], dtype=np.int64)
        shared[-1] = interface
        strip_ids = attach(strip, shared)
        blocks.append((3, quads_of(strip_ids, False)))
        gfaces["ground"] = np.stack([strip_ids[0, :-1], strip_ids[0, 1:]], axis=1)
        for name, column in (("inlet", 0), ("outlet", -1)):
            edges = np.stack([strip_ids[:-1, column], strip_ids[1:, column]], axis=1)
            gfaces[name] = np.concatenate([gfaces[name], edges])
    points = np.concatenate(points)
    kinds = {"ground": "wall" if ground_wall else "symmetry", "inlet": "velocity-inlet",
             "outlet": "pressure-outlet", "top": "symmetry"}
    boundaries += [(k, kinds[k], v) for k, v in gfaces.items()]
    if png:
        preview(png, points, blocks, [(-0.3, 1.5, ground_y - 0.05, 0.4), (te[0] - 0.02, te[0] + 0.04, te[1] - 0.025, te[1] + 0.025),
                                      (-0.03, 0.05, -0.05, 0.03)])
    stats = fluent_msh.write(out, points, blocks, boundaries, comment=f"fsae_aero wing2d {meta['case']}")
    counts = {}
    for t, b in blocks:
        key = "tri" if t == 2 else "quad"
        counts[key] = counts.get(key, 0) + int(len(b))
    meta.update(stats, level=level, blunt_te=bool(blunt), ds_wall=ds_wall, growth=ratio, oblock_layers=nl - 1,
                oblock_thickness=float(np.max(np.linalg.norm(block[-1] - block[0], axis=1))), wing_faces=int(nloop),
                ground_y=ground_y, gap=gap, cell_types=counts,
                domain={"x_in": x_in, "x_out": x_out, "y_far": y_far})
    if len(ground_offsets) > 1:
        meta.update(ground_layers=len(ground_offsets) - 1, ground_first_height=ds_wall,
                    ground_layer_thickness=ground_offsets[-1], ground_far_layers="uniform, thickness 1c beyond [-4, 7]")
    if wake_refinement != 1:
        meta["wake_refinement"] = wake_refinement
        meta["wake_refinement_half_width"] = 0.45
        meta["wake_refinement_length"] = 3.0
    if trailing_edge_refinement != 1:
        meta["trailing_edge_refinement"] = trailing_edge_refinement
        meta["trailing_edge_refinement_start"] = 0.85
        meta["trailing_edge_refinement_coordinate"] = "chord projection"
    if wall_refinement != 1:
        meta["wall_refinement"] = wall_refinement
    pathlib.Path(str(out) + ".json").write_text(json.dumps(meta, indent=1))
    return meta


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case", choices=["zerihan", "s1223"])
    ap.add_argument("--hc", default="free", help="ride height h/c, or 'free'")
    ap.add_argument("--alpha", type=float, default=0.0, help="S1223 angle of attack, deg")
    ap.add_argument("--level", type=float, default=2, help="fractional levels interpolate the refinement factor")
    ap.add_argument("--out", required=True)
    ap.add_argument("--png", help="write a mesh preview image")
    ap.add_argument("--far", type=float, default=1.0, help="scale the far-field distances")
    ap.add_argument("--ground-layers", action="store_true", help="resolve the ground with a conformal quad strip")
    ap.add_argument("--domain", type=float, nargs=3, metavar=("X_IN", "X_OUT", "Y_TOP"),
                    help="explicit domain boundaries in chords (overrides --far)")
    ap.add_argument("--wake-refinement", type=float, default=1.0, help="divide wake-box cell sizes by this factor")
    ap.add_argument("--trailing-edge-refinement", type=int, default=1, help="split wall faces aft of x/c .85 and the base")
    ap.add_argument("--wall-refinement", type=float, default=1.0, help="refine surface/wall layers while keeping outer sizes fixed")
    a = ap.parse_args()
    if not math.isfinite(a.level) or a.level <= 0:
        ap.error("--level must be finite and positive")
    if a.ground_layers and (a.case != "zerihan" or a.hc == "free"):
        ap.error("--ground-layers requires zerihan with a numeric ride height")
    if a.domain and not (all(math.isfinite(v) for v in a.domain) and
                         a.domain[0] < -1 and a.domain[1] > 2 and a.domain[2] > 1):
        ap.error("--domain requires finite X_IN < -1, X_OUT > 2 and Y_TOP > 1")
    pivot = np.array([0.25, 0.0])
    if a.case == "zerihan":
        contour = zerihan_contour()
        pivot = np.array([0.25 * contour[-1][0], 0.25 * 0.5 * (contour[0][1] + contour[-1][1])])
        meta = {"case": "zerihan", "alpha_ref_deg": 1.0, "true_incidence_deg": 3.45, "hc": a.hc}
        if a.hc == "free":
            ground_y, wall = -Y_FAR * a.far, False
        else:
            # h: ground to the lowest suction-surface point with the wing at the
            # zero reference incidence (rotate -1 deg about c/4), SOURCE.md.
            low = rotate(contour, -1.0, pivot)[:, 1].min()
            ground_y, wall = low - float(a.hc), True
            meta["lowest_point_y_alpha0"] = float(low)
    else:
        contour = selig_contour()
        contour = rotate(contour, -a.alpha, pivot)
        meta = {"case": "s1223", "alpha_deg": a.alpha, "coordinates": "as-tested (dap)"}
        ground_y, wall = -Y_FAR * a.far, False
    pathlib.Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    if not math.isfinite(a.wake_refinement) or a.wake_refinement <= 0:
        ap.error("--wake-refinement must be finite and positive")
    if a.trailing_edge_refinement < 1:
        ap.error("--trailing-edge-refinement must be at least one")
    if not math.isfinite(a.wall_refinement) or a.wall_refinement <= 0:
        ap.error("--wall-refinement must be finite and positive")
    print(json.dumps(build(contour, ground_y, a.level, a.out, wall, meta, a.png, a.far, ground_layers=a.ground_layers,
                           domain=a.domain, wake_refinement=a.wake_refinement,
                           trailing_edge_refinement=a.trailing_edge_refinement, wall_refinement=a.wall_refinement)))


if __name__ == "__main__":
    main()
