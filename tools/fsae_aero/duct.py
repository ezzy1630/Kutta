#!/usr/bin/env python3
"""Minimal reproduction: the solver diverges on gmsh tetrahedra (REPORT.md).

A 2 x 1 x 1 duct (velocity inlet at x = 0, pressure outlet at x = 2, walls
elsewhere) meshed three ways at the same nominal size:

  tet   gmsh Delaunay tetrahedra
  hex   gmsh transfinite hexahedra
  split the tetrahedra split into hexahedra by gmsh (SubdivisionAlgorithm 2)

  build-py/bin/python tools/fsae_aero/duct.py tet --out /tmp/duct_tet.msh
  build/bench_fsae_aero --mesh /tmp/duct_tet.msh --re 1e5 --walls walls \
      --backend cpu --iterations 200 --check 20
"""
import argparse
import json
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import fluent_msh  # noqa: E402

NODES = {2: 3, 3: 4, 4: 4, 5: 8, 6: 6, 7: 5}  # gmsh element type -> node count


def build(kind, size, out):
    import gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    occ = gmsh.model.occ
    box = occ.addBox(0, 0, 0, 2, 1, 1)
    occ.synchronize()
    groups = {"inlet": [], "outlet": [], "walls": []}
    for _, s in gmsh.model.getBoundary([(3, box)], oriented=False):
        x0, _, _, x1, _, _ = gmsh.model.getBoundingBox(2, s)
        groups["inlet" if x1 < 1e-6 else "outlet" if x0 > 2 - 1e-6 else "walls"].append(s)
    gmsh.option.setNumber("Mesh.MeshSizeMin", size)
    gmsh.option.setNumber("Mesh.MeshSizeMax", size)
    if kind == "hex":
        for _, c in gmsh.model.getEntities(1):
            gmsh.model.mesh.setTransfiniteCurve(c, int(round(occ.getMass(1, c) / size)) + 1)
        for _, s in gmsh.model.getEntities(2):
            gmsh.model.mesh.setTransfiniteSurface(s)
            gmsh.model.mesh.setRecombine(2, s)
        gmsh.model.mesh.setTransfiniteVolume(box)
    if kind == "split":
        gmsh.option.setNumber("Mesh.SubdivisionAlgorithm", 2)
    gmsh.model.mesh.generate(3)
    tags, coords, _ = gmsh.model.mesh.getNodes()
    points = np.zeros((int(tags.max()) + 1, 3))
    points[tags.astype(np.int64)] = coords.reshape(-1, 3)
    types, _, nodes = gmsh.model.mesh.getElements(3)
    blocks = [(int(t), n.astype(np.int64).reshape(-1, NODES[int(t)])) for t, n in zip(types, nodes)]
    kinds = {"inlet": "velocity-inlet", "outlet": "pressure-outlet", "walls": "wall"}
    boundary = []
    for name, surfaces in groups.items():
        faces = []
        for s in surfaces:
            types, _, nodes = gmsh.model.mesh.getElements(2, s)
            for t, n in zip(types, nodes):
                f = n.astype(np.int64).reshape(-1, NODES[int(t)])
                faces.append(np.concatenate([f, -np.ones((len(f), 4 - f.shape[1]), dtype=np.int64)], axis=1))
        boundary.append((name, kinds[name], np.concatenate(faces)))
    gmsh.finalize()
    return fluent_msh.write(out, points, blocks, boundary, comment=f"fsae_aero duct {kind}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("kind", choices=["tet", "hex", "split"])
    ap.add_argument("--size", type=float, default=0.1)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    print(json.dumps(build(a.kind, a.size, a.out)))


if __name__ == "__main__":
    main()
