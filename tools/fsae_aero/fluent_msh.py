"""Write a conformal gmsh mesh as an ASCII Fluent mesh (.msh) that Kutta reads.

Only what Kutta's reader needs is written: nodes, one fluid cell zone, an
interior face zone and one face zone per named boundary. Kutta's conventions
(src/kutta/io/fluent_common.cpp): in 3D a face's node winding points into c0;
in 2D the directed edge n0->n1 has c0 on its left. Orientation is derived
geometrically, so gmsh's local node orders never decide a sign.
"""
import numpy as np

# gmsh element type -> (dimension, Fluent element code, local faces)
CELLS = {
    2: (2, 1, [(0, 1), (1, 2), (2, 0)]),
    3: (2, 3, [(0, 1), (1, 2), (2, 3), (3, 0)]),
    4: (3, 2, [(0, 2, 1), (0, 1, 3), (0, 3, 2), (1, 2, 3)]),
    5: (3, 4, [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]),
    6: (3, 6, [(0, 2, 1), (3, 4, 5), (0, 1, 4, 3), (1, 2, 5, 4), (2, 0, 3, 5)]),
    7: (3, 5, [(0, 3, 2, 1), (0, 1, 4), (1, 2, 4), (2, 3, 4), (3, 0, 4)]),
}
BOUNDARY_TYPES = {"wall": 3, "pressure-outlet": 5, "symmetry": 7, "velocity-inlet": 10}


def _cell_centroids(points, blocks):
    return [points[nodes].mean(axis=1) for _, nodes in blocks]


def _faces(points, blocks, dim):
    """All cell faces as (nodes padded with -1, cell id); outward for 3D, cell-left for 2D."""
    out = []
    first = 0
    orientation = 0
    for (etype, nodes), centre in zip(blocks, _cell_centroids(points, blocks)):
        _, _, local = CELLS[etype]
        cells = np.arange(first, first + len(nodes))
        if dim == 2:
            # Cells must already share one orientation; flipping an inverted
            # cell would hide a tangled mesh. Counter-clockwise keeps the cell
            # on each edge's left.
            x, y = points[nodes, 0], points[nodes, 1]
            area = 0.5 * np.sum(x * np.roll(y, -1, axis=1) - np.roll(x, -1, axis=1) * y, axis=1)
            sign = orientation if orientation else np.sign(np.sum(np.sign(area)))
            orientation = sign
            if np.any(area * sign <= 0):
                bad = np.flatnonzero(area * sign <= 0)
                where = np.round(points[nodes[bad[:5]]].mean(axis=1)[:, :2], 6).tolist()
                raise ValueError(f"{len(bad)} inverted or degenerate 2D cells, e.g. at {where}")
            if sign < 0:
                nodes = nodes[:, ::-1]
            k = nodes.shape[1]
            for a in range(k):
                out.append((np.stack([nodes[:, a], nodes[:, (a + 1) % k]], axis=1), cells))
        else:
            for face in local:
                f = nodes[:, face]
                p = points[f]
                # Newell normal: exact for planar polygons, robust for warped quads.
                normal = np.zeros((len(f), 3))
                for a in range(len(face)):
                    normal += np.cross(p[:, a], p[:, (a + 1) % len(face)])
                outward = np.einsum("ij,ij->i", normal, p.mean(axis=1) - centre) > 0
                f = np.where(outward[:, None], f, f[:, ::-1])
                if len(face) == 3:
                    f = np.concatenate([f, -np.ones((len(f), 1), dtype=f.dtype)], axis=1)
                out.append((f, cells))
        first += len(nodes)
    nodes = np.concatenate([f for f, _ in out])
    cells = np.concatenate([c for _, c in out])
    return nodes, cells


def _key(nodes):
    k = np.sort(nodes, axis=1)
    return k


def write(path, points, blocks, boundaries, comment="kutta fsae_aero"):
    """points: (N,3) array; blocks: [(gmsh type, (n,k) 0-based node array)];
    boundaries: [(name, fluent type name, (m,k) node array)] in zone order."""
    dim = CELLS[blocks[0][0]][0]
    if any(CELLS[t][0] != dim for t, _ in blocks):
        raise ValueError("mixed-dimension cell blocks")
    used = np.unique(np.concatenate([b.ravel() for _, b in blocks]))
    remap = -np.ones(len(points), dtype=np.int64)
    remap[used] = np.arange(len(used))
    points = points[used]
    blocks = [(t, remap[b]) for t, b in blocks]
    original = boundaries
    # -1 pads triangles among quads; keep it out of the remap.
    boundaries = [(n, t, np.where(b < 0, -1, remap[b])) for n, t, b in boundaries]
    if any(np.any((b < 0) & (b0 >= 0)) for (_, _, b), (_, _, b0) in zip(boundaries, original)):
        raise ValueError("boundary face references a node no cell uses")

    nodes, cells = _faces(points, blocks, dim)
    key = _key(nodes)
    order = np.lexsort(key.T[::-1])
    key, nodes, cells = key[order], nodes[order], cells[order]
    same = np.all(key[1:] == key[:-1], axis=1)
    if np.any(same[1:] & same[:-1]):
        raise ValueError("a face is shared by more than two cells")
    first = np.ones(len(key), dtype=bool)
    first[1:] = ~same
    pair = np.zeros(len(key), dtype=bool)
    pair[:-1] = same
    starts = np.flatnonzero(first)
    interior = starts[pair[starts]]
    single = starts[~pair[starts]]
    # 3D: nodes are outward from the first cell, so they point into the second (c0).
    # 2D: the first cell is on the edge's left, so it is c0.
    if dim == 3:
        inner = (nodes[interior], cells[interior + 1], cells[interior])
    else:
        inner = (nodes[interior], cells[interior], cells[interior + 1])

    lookup = {tuple(k): i for i, k in enumerate(key[single])}
    zone_rows = []
    claimed = np.zeros(len(single), dtype=bool)
    for name, kind, faces in boundaries:
        padded = faces if faces.shape[1] == key.shape[1] else np.concatenate(
            [faces, -np.ones((len(faces), key.shape[1] - faces.shape[1]), dtype=faces.dtype)], axis=1)
        rows = np.array([lookup[tuple(k)] for k in _key(padded)], dtype=np.int64)
        if np.any(claimed[rows]):
            raise ValueError("boundary zones overlap: " + name)
        claimed[rows] = True
        f = nodes[single[rows]]
        c = cells[single[rows]]
        # 3D boundary winding must point into the only cell; 2D keeps the cell on the left.
        if dim == 3:
            f = np.where((f[:, 3] < 0)[:, None], f[:, [2, 1, 0, 3]], f[:, ::-1])
        zone_rows.append((name, kind, f, c))
    if not np.all(claimed):
        raise ValueError(f"{np.count_nonzero(~claimed)} boundary faces have no zone")

    n_faces = len(interior) + len(single)
    n_cells = sum(len(b) for _, b in blocks)
    with open(path, "w") as out:
        out.write(f'(0 "{comment}")\n(2 {dim})\n')
        out.write(f"(10 (0 1 {len(points):x} 0 {dim}))\n(10 (1 1 {len(points):x} 1 {dim})(\n")
        np.savetxt(out, points[:, :dim], fmt="%.17g")
        out.write("))\n")
        out.write(f"(12 (0 1 {n_cells:x} 0))\n(12 (2 1 {n_cells:x} 1 0)(\n")
        np.savetxt(out, np.concatenate([np.full(len(b), CELLS[t][1]) for t, b in blocks]), fmt="%x")
        out.write("))\n")
        out.write(f"(13 (0 1 {n_faces:x} 0))\n")
        start = 1
        zones = [(3, "interior", "default-interior", inner)]
        for z, (name, kind, f, c) in enumerate(zone_rows):
            zones.append((4 + z, kind, name, (f, c, np.full(len(c), -1))))
        for zid, kind, name, (f, c0, c1) in zones:
            code = 2 if kind == "interior" else BOUNDARY_TYPES[kind]
            end = start + len(f) - 1
            out.write(f"(13 ({zid:x} {start:x} {end:x} {code:x} 0)(\n")
            # Face ids within a zone are arbitrary, so group rows by node count.
            count = np.count_nonzero(f >= 0, axis=1)
            for k in np.unique(count):
                sel = count == k
                rows = np.concatenate([np.full((sel.sum(), 1), k), f[sel][:, :k] + 1,
                                       c0[sel][:, None] + 1, c1[sel][:, None] + 1], axis=1)
                np.savetxt(out, rows, fmt="%x")
            out.write("))\n")
            start = end + 1
        out.write("(45 (2 fluid fluid)())\n")
        for zid, kind, name, _ in zones:
            out.write(f"(45 ({zid} {kind} {name})())\n")  # decimal, unlike section headers
    return {"cells": n_cells, "faces": n_faces, "nodes": len(points)}
