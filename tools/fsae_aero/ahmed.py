#!/usr/bin/env python3
"""Ahmed body meshes (Fluent .msh) for the FSAE aero validation.

Geometry from the documented dimensions (ERCOFTAC AC1-05; Ahmed et al. 1984):
L = 1.044 m, W = 0.389 m, H = 0.288 m, 100 mm radii on the four front edges,
222 mm slant at `--slant` degrees, sharp rear edges, 50 mm ground clearance.
Stilts are omitted, as in the CFX-5 SST reference (CFD11). Coordinates follow
Lienhart & Becker: x = 0 at the base, y = 0 mid-plane, z = 0 ground. The full
body is meshed (no symmetry plane) in CFD11's domain mirrored about y = 0:
[-3L, 6L] from the nose x [-0.9L, 0.9L] x [0, 1.15L]. The ground is slip ahead
of 2L before the nose and a stationary wall after it.

Mesh: gmsh surface triangles on the body, prism layers marched outward from
them (first height `FIRST`), gmsh tetrahedra outside the
layer shell, then every cell split into hexahedra (subdivide.py), because the
solver diverged on tetrahedra when this was written (fsae_aero/REPORT.md).

Zones: body (wall), ground (wall), ground_slip, side, top (symmetry),
inlet (velocity), outlet (pressure).
"""
import argparse
import json
import math
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import fluent_msh  # noqa: E402
import subdivide  # noqa: E402

L, W, H, CLEAR, R, SLANT = 1.044, 0.389, 0.288, 0.050, 0.100, 0.222
# 50 um is y+ of a few at 40 m/s. The CFX-5 reference resolves y+ ~ 1
# (--first 1.5e-5). Thin layers diverged at the sharp rear edges until the
# poor-cell pressure treatment (external_3d_robustness/REPORT.md); the default
# stays 50 um so the existing meshes and tests keep their names and values.
FIRST, GROWTH = 5e-5, 1.35


def body(occ, slant_deg):
    import gmsh
    phi = math.radians(slant_deg)
    b = occ.addBox(-L, -W / 2, CLEAR, L, W, H)
    occ.synchronize()
    front = [e for _, e in gmsh.model.getEntitiesInBoundingBox(
        -L - 1e-6, -W, CLEAR - 1e-6, -L + 1e-6, W, CLEAR + H + 1e-6, 1)]
    b = occ.fillet([b], front, [R], removeVolume=True)[0][1]
    xs, top = -SLANT * math.cos(phi), CLEAR + H
    cutter = occ.addBox(xs, -W, top, 1.0, 2 * W, 1.0)
    occ.rotate([(3, cutter)], xs, 0, top, 0, 1, 0, phi)
    out, _ = occ.cut([(3, b)], [(3, cutter)])
    return out[0][1]


def surface(slant, size):
    """Outward-oriented body surface triangles (points, triangles)."""
    import gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("body")
    vb = body(gmsh.model.occ, slant)
    gmsh.model.occ.synchronize()
    for key, value in {"Mesh.MeshSizeMax": size, "Mesh.MeshSizeMin": 0.25 * size,
                       "Mesh.MeshSizeFromCurvature": 30, "Mesh.Algorithm": 6}.items():
        gmsh.option.setNumber(key, value)
    gmsh.model.mesh.generate(2)
    tags, coords, _ = gmsh.model.mesh.getNodes()
    points = np.zeros((int(tags.max()) + 1, 3))
    points[tags.astype(np.int64)] = coords.reshape(-1, 3)
    tris = np.concatenate([gmsh.model.mesh.getElements(2, s)[2][0].reshape(-1, 3)
                           for _, s in gmsh.model.getBoundary([(3, vb)], oriented=False)]).astype(np.int64)
    gmsh.finalize()
    used = np.unique(tris)
    remap = -np.ones(len(points), dtype=np.int64)
    remap[used] = np.arange(len(used))
    points, tris = points[used], remap[tris]
    # The body is convex, so outward means away from an interior point.
    n = np.cross(points[tris[:, 1]] - points[tris[:, 0]], points[tris[:, 2]] - points[tris[:, 0]])
    inside = np.array([-L / 2, 0, CLEAR + H / 2])
    flip = np.einsum("ij,ij->i", n, points[tris].mean(axis=1) - inside) < 0
    tris[flip] = tris[flip][:, ::-1]
    return points, tris


def vertex_normals(points, tris):
    n = np.cross(points[tris[:, 1]] - points[tris[:, 0]], points[tris[:, 2]] - points[tris[:, 0]])
    acc = np.zeros_like(points)
    for k in range(3):
        np.add.at(acc, tris[:, k], n)
    return acc / np.linalg.norm(acc, axis=1)[:, None]


def layers(points, tris, top_size, first):
    """March prism layers: geometric heights from `first`, capped at `top_size`
    and stopped below 40% of the ground clearance."""
    shells, h, total = [points], first, 0.0
    while total + h <= 0.4 * CLEAR and h <= top_size:
        shells.append(shells[-1] + h * vertex_normals(shells[-1], tris))
        total += h
        h *= GROWTH
    return np.array(shells), total


def build(slant, level, out, first=FIRST):
    import gmsh
    r = 2.0 ** ((level - 2) / 2)  # surface spacing /sqrt(2), cells ~x2 per level
    size = 0.028 / r
    points, tris = surface(slant, size)
    shells, thickness = layers(points, tris, 0.6 * size, first)
    nl, nn = shells.shape[0], shells.shape[1]
    top = shells[-1]

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("ahmed")
    x0, x1, ymax, zmax, xg = -L - 3 * L, -L + 6 * L, 0.9 * L, 1.15 * L, -3 * L
    geo = gmsh.model.geo
    corner = {}
    for x in (x0, xg, x1):
        for y in (-ymax, ymax):
            for z in (0.0, zmax):
                if x == xg and z > 0:
                    continue
                corner[x, y, z] = geo.addPoint(x, y, z)
    line = {}

    def ln(a, b):
        key = (a, b)
        if key not in line:
            line[key] = geo.addLine(corner[a], corner[b])
            line[b, a] = -line[key]
        return line[key]

    def face(*pts):
        loop = geo.addCurveLoop([ln(pts[i], pts[(i + 1) % len(pts)]) for i in range(len(pts))])
        return geo.addPlaneSurface([loop])

    lo, hi = -ymax, ymax
    groups = {
        "inlet": [face((x0, lo, 0.0), (x0, hi, 0.0), (x0, hi, zmax), (x0, lo, zmax))],
        "outlet": [face((x1, lo, 0.0), (x1, lo, zmax), (x1, hi, zmax), (x1, hi, 0.0))],
        "ground_slip": [face((x0, lo, 0.0), (xg, lo, 0.0), (xg, hi, 0.0), (x0, hi, 0.0))],
        "ground": [face((xg, lo, 0.0), (x1, lo, 0.0), (x1, hi, 0.0), (xg, hi, 0.0))],
        "top": [face((x0, lo, zmax), (x0, hi, zmax), (x1, hi, zmax), (x1, lo, zmax))],
        "side": [face((x0, lo, 0.0), (x0, lo, zmax), (x1, lo, zmax), (x1, lo, 0.0), (xg, lo, 0.0)),
                 face((x0, hi, 0.0), (xg, hi, 0.0), (x1, hi, 0.0), (x1, hi, zmax), (x0, hi, zmax))],
    }
    geo.synchronize()
    # The layer shell's outer surface, as an already-meshed discrete surface.
    shell = gmsh.model.addDiscreteEntity(2)
    first_tag = 10_000_000
    gmsh.model.mesh.addNodes(2, shell, first_tag + np.arange(nn), top.ravel())
    gmsh.model.mesh.addElementsByType(shell, 2, [], (first_tag + tris).ravel())
    outer = geo.addSurfaceLoop([s for v in groups.values() for s in v])
    inner = geo.addSurfaceLoop([shell])
    vol = geo.addVolume([outer, inner])
    geo.synchronize()
    f = gmsh.model.mesh.field
    dist = f.add("Distance")
    f.setNumbers(dist, "SurfacesList", [shell])
    thr = f.add("Threshold")
    f.setNumber(thr, "InField", dist)
    f.setNumber(thr, "SizeMin", size)
    f.setNumber(thr, "SizeMax", 0.4)
    f.setNumber(thr, "DistMin", 0.0)
    f.setNumber(thr, "DistMax", 3.0)
    wake = f.add("Box")
    f.setNumber(wake, "VIn", 1.6 * size)
    f.setNumber(wake, "VOut", 0.4)
    for k, v in {"XMin": -L - 0.1, "XMax": 1.0, "YMin": -0.4, "YMax": 0.4, "ZMin": 0, "ZMax": 0.5, "Thickness": 2.0}.items():
        f.setNumber(wake, k, v)
    fmin = f.add("Min")
    f.setNumbers(fmin, "FieldsList", [thr, wake])
    f.setAsBackgroundMesh(fmin)
    for key, value in {"Mesh.MeshSizeExtendFromBoundary": 0, "Mesh.MeshSizeFromPoints": 0,
                       "Mesh.MeshSizeFromCurvature": 0, "Mesh.Algorithm": 6, "Mesh.Algorithm3D": 1,
                       "Mesh.Optimize": 1, "General.NumThreads": 8}.items():
        gmsh.option.setNumber(key, value)
    gmsh.model.mesh.generate(3)
    tags, coords, _ = gmsh.model.mesh.getNodes()
    gpts = np.zeros((int(tags.max()) + 1, 3))
    gpts[tags.astype(np.int64)] = coords.reshape(-1, 3)
    ts, _, ns = gmsh.model.mesh.getElements(3, vol)
    tets = np.concatenate([n.astype(np.int64).reshape(-1, 4) for t, n in zip(ts, ns) if int(t) == 4])
    if sum(len(n) for n in ns) != 4 * len(tets):
        raise SystemExit("gmsh produced non-tetrahedral volume elements")
    bnd = {}
    for name, surfaces in groups.items():
        bnd[name] = np.concatenate([gmsh.model.mesh.getElements(2, s)[2][0].astype(np.int64).reshape(-1, 3)
                                    for s in surfaces])
    # gmsh renumbers nodes when meshing; find the shell nodes by their exact coordinates.
    stags, scoords, _ = gmsh.model.mesh.getNodes(2, shell, includeBoundary=True)
    gmsh.finalize()
    where = {tuple(c): int(t) for t, c in zip(stags, scoords.reshape(-1, 3).tolist())}
    try:
        shell_ids = np.array([where[tuple(c)] for c in top.tolist()], dtype=np.int64)
    except KeyError:
        raise SystemExit("gmsh moved layer-shell nodes")

    # Layers 0..nl-2 get new ids after gmsh's; the outermost shell keeps gmsh's.
    ids = np.empty((nl, nn), dtype=np.int64)
    ids[:-1] = len(gpts) + np.arange((nl - 1) * nn).reshape(nl - 1, nn)
    ids[-1] = shell_ids
    allpts = np.concatenate([gpts, shells[:-1].reshape(-1, 3)])
    prisms = np.concatenate([np.concatenate([ids[k][tris], ids[k + 1][tris]], axis=1) for k in range(nl - 1)])
    kinds = {"body": "wall", "ground": "wall", "ground_slip": "symmetry", "side": "symmetry", "top": "symmetry",
             "inlet": "velocity-inlet", "outlet": "pressure-outlet"}
    pad = lambda f: np.concatenate([f, -np.ones((len(f), 1), dtype=np.int64)], axis=1)
    boundary = [("body", "wall", pad(ids[0][tris]))] + [(k, kinds[k], pad(v)) for k, v in bnd.items()]
    hpts, hexes, hbnd = subdivide.subdivide(allpts, tets, prisms, boundary)
    stats = fluent_msh.write(out, hpts, [(5, hexes)], hbnd, comment=f"fsae_aero ahmed {slant}")
    meta = dict(stats, case="ahmed", slant_deg=slant, level=level, surface_size=size, body_triangles=int(len(tris)),
                prism_layers=int(nl - 1), layer_thickness=thickness, first_height=first, growth=GROWTH,
                tets=int(len(tets)), prisms=int(len(prisms)), hexes=int(len(hexes)),
                domain={"x": [x0, x1], "y": [-ymax, ymax], "z": [0, zmax], "ground_wall_from_x": xg})
    pathlib.Path(str(out) + ".json").write_text(json.dumps(meta, indent=1))
    return meta


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--slant", type=float, default=35.0)
    ap.add_argument("--level", type=int, default=2)
    ap.add_argument("--out", required=True)
    ap.add_argument("--first", type=float, default=FIRST, help="first prism height, m")
    a = ap.parse_args()
    pathlib.Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    print(json.dumps(build(a.slant, a.level, a.out, a.first)))


if __name__ == "__main__":
    main()
