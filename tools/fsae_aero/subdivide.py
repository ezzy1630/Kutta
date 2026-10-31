"""Conformal all-hex subdivision of a tet + prism mesh.

Kutta's current solver diverges on gmsh tetrahedra (reported in
validation/results/fsae_aero/REPORT.md) but is stable on hexahedra, so tets are
split into four hexes (vertex, edge midpoints, face centroids, cell centroid).
Boundary-layer prisms are split into three hexes through their full height,
which keeps the wall-normal spacing (y+) of the layer unchanged. Shared edges
and triangles get one midpoint/centroid each, so the result stays conformal.
"""
import numpy as np


class _Points:
    def __init__(self, points):
        self.chunks = [points]
        self.count = len(points)

    def add(self, pts):
        ids = self.count + np.arange(len(pts))
        self.chunks.append(pts)
        self.count += len(pts)
        return ids

    def array(self):
        return np.concatenate(self.chunks)


def _unique_ids(keys, values, points):
    """Map each row of `keys` (sorted node tuples) to a new point (mean of `values`)."""
    uniq, inverse = np.unique(keys, axis=0, return_inverse=True)
    first = np.zeros(len(uniq), dtype=np.int64)
    first[inverse[::-1]] = np.arange(len(keys))[::-1]
    return uniq, inverse.ravel(), points.add(values[first])


def subdivide(points, tets, prisms, boundary):
    """points (N,3); tets (n,4); prisms (m,6) with node i+3 above node i;
    boundary: [(name, kind, faces padded to 4 with -1)]. Returns hex mesh."""
    p = points.astype(float)
    store = _Points(p)
    # Edges needing midpoints: all tet edges and prism triangle edges.
    tri_e = [(0, 1), (1, 2), (2, 0)]
    tet_e = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    edges = [tets[:, list(e)] for e in tet_e] + [prisms[:, [a + o, b + o]] for o in (0, 3) for a, b in tri_e]
    edges = np.sort(np.concatenate(edges), axis=1)
    e_uniq, _, e_ids = _unique_ids(edges, p[edges].mean(axis=1), store)
    e_lookup = {tuple(k): i for k, i in zip(e_uniq.tolist(), e_ids.tolist())}
    # Triangles needing centroids: all tet faces and prism top/bottom.
    tet_f = [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]
    faces = [tets[:, list(f)] for f in tet_f] + [prisms[:, [0, 1, 2]], prisms[:, [3, 4, 5]]]
    faces = np.sort(np.concatenate(faces), axis=1)
    f_uniq, _, f_ids = _unique_ids(faces, p[faces].mean(axis=1), store)
    f_lookup = {tuple(k): i for k, i in zip(f_uniq.tolist(), f_ids.tolist())}

    def mid(a, b):
        k = np.sort(np.stack([a, b], axis=1), axis=1)
        return np.array([e_lookup[t] for t in map(tuple, k.tolist())], dtype=np.int64)

    def cen(a, b, c):
        k = np.sort(np.stack([a, b, c], axis=1), axis=1)
        return np.array([f_lookup[t] for t in map(tuple, k.tolist())], dtype=np.int64)

    hexes = []
    if len(tets):
        centre = store.add(p[tets].mean(axis=1))
        for v in range(4):
            a, b, c = [w for w in range(4) if w != v]
            V, A, B, C = (tets[:, i] for i in (v, a, b, c))
            hexes.append(np.stack([V, mid(V, A), cen(V, A, B), mid(V, B),
                                   mid(V, C), cen(V, A, C), centre, cen(V, B, C)], axis=1))
    if len(prisms):
        lo, hi = cen(*(prisms[:, i] for i in (0, 1, 2))), cen(*(prisms[:, i] for i in (3, 4, 5)))
        for v in range(3):
            a, b = (v + 1) % 3, (v + 2) % 3
            V, A, B = (prisms[:, i] for i in (v, a, b))
            V2, A2, B2 = (prisms[:, i + 3] for i in (v, a, b))
            hexes.append(np.stack([V, mid(V, A), lo, mid(V, B), V2, mid(V2, A2), hi, mid(V2, B2)], axis=1))
    hexes = np.concatenate(hexes)

    tri_edges = set()
    for o in (0, 3):
        for a, b in tri_e:
            tri_edges.update(map(tuple, np.sort(prisms[:, [a + o, b + o]], axis=1).tolist()))
    out = []
    for name, kind, f in boundary:
        tri, quad = f[f[:, 3] < 0][:, :3], f[f[:, 3] >= 0]
        pieces = []
        if len(tri):
            m = [mid(tri[:, i], tri[:, (i + 1) % 3]) for i in range(3)]
            c = cen(tri[:, 0], tri[:, 1], tri[:, 2])
            for i in range(3):
                pieces.append(np.stack([tri[:, i], m[i], c, m[(i + 2) % 3]], axis=1))
        if len(quad):
            # Prism side faces: split through the midpoints of the two triangle edges.
            first = np.array([tuple(sorted(e)) in tri_edges for e in quad[:, :2].tolist()])
            q = np.where(first[:, None], quad, np.roll(quad, -1, axis=1))
            m01, m23 = mid(q[:, 0], q[:, 1]), mid(q[:, 2], q[:, 3])
            pieces += [np.stack([q[:, 0], m01, m23, q[:, 3]], axis=1), np.stack([m01, q[:, 1], q[:, 2], m23], axis=1)]
        out.append((name, kind, np.concatenate(pieces)))
    return store.array(), hexes, out
