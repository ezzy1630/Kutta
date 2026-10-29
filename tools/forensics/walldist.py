"""Own wall-distance candidates, compared against Fluent's saved SV_WALL_DIST.

d_fc    : distance from cell centroid to nearest wall-face centroid
d_exact : exact Euclidean distance from cell centroid to the nearest wall polygon
          (each polygon fanned into triangles about its node average)
d_plane : |(x - x_f) . n_f| to the plane of the face that is nearest in d_exact
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.spatial import cKDTree

import cff


def wall_faces(mesh):
    return np.concatenate([np.arange(z["lo"], z["hi"]) for z in mesh["zones"] if z["type"] == 3])


def point_triangle_distance(p, a, b, c):
    """Vectorised exact distance from points p to triangles (a, b, c) (Ericson, RTCD 5.1.5)."""
    ab, ac, ap = b - a, c - a, p - a
    d1, d2 = (ab * ap).sum(-1), (ac * ap).sum(-1)
    bp = p - b
    d3, d4 = (ab * bp).sum(-1), (ac * bp).sum(-1)
    cp = p - c
    d5, d6 = (ab * cp).sum(-1), (ac * cp).sum(-1)
    va = d3 * d6 - d5 * d4
    vb = d5 * d2 - d1 * d6
    vc = d1 * d4 - d3 * d2
    denom = va + vb + vc
    with np.errstate(divide="ignore", invalid="ignore"):
        v = vb / denom
        w = vc / denom
        q = a + ab * v[..., None] + ac * w[..., None]  # interior
        t_ab = np.clip(d1 / (d1 - d3), 0, 1)
        t_ac = np.clip(d2 / (d2 - d6), 0, 1)
        t_bc = np.clip((d4 - d3) / ((d4 - d3) + (d5 - d6)), 0, 1)
    cand = [q, a, b, c, a + ab * t_ab[..., None], a + ac * t_ac[..., None], b + (c - b) * t_bc[..., None]]
    # Taking the min over all candidates is safe: every candidate is a point on the triangle,
    # and the true closest point is one of them (interior only valid when inside).
    inside = (va >= 0) & (vb >= 0) & (vc >= 0)
    dist = np.full(p.shape[:-1], np.inf)
    for i, qq in enumerate(cand):
        dd = np.linalg.norm(p - qq, axis=-1)
        if i == 0:
            dd = np.where(inside, dd, np.inf)
        if i >= 4:
            dd = np.where(np.isfinite(dd), dd, np.inf)
        dist = np.minimum(dist, dd)
    return dist


def wall_distances(mesh, k_candidates=24):
    wf = wall_faces(mesh)
    nodes, nn, fn = mesh["nodes"], mesh["nnodes"], mesh["fnodes"]
    ptr = np.concatenate([[0], np.cumsum(nn)])
    fc = mesh["fcent"][wf]
    nhat = mesh["farea"][wf] / np.linalg.norm(mesh["farea"][wf], axis=1)[:, None]
    cen = mesh["cen"]
    tree = cKDTree(fc)
    d_fc, idx = tree.query(cen, k=k_candidates)
    best = np.full(cen.shape[0], np.inf)
    best_face = np.zeros(cen.shape[0], dtype=np.int64)
    # Triangle fan per wall face: (node_i, node_{i+1}, node average).
    maxn = nn[wf].max()
    for j in range(k_candidates):
        f_local = idx[:, j]
        f = wf[f_local]
        avg = np.zeros((f.size, 3))
        for s in range(maxn):
            has = s < nn[f]
            avg[has] += nodes[fn[ptr[f[has]] + s]]
        avg /= nn[f][:, None]
        dmin = np.full(f.size, np.inf)
        for s in range(maxn):
            has = s < nn[f]
            if not has.any():
                continue
            a = nodes[fn[ptr[f[has]] + s]]
            b = nodes[fn[ptr[f[has]] + (s + 1) % nn[f[has]]]]
            dd = point_triangle_distance(cen[has], a, b, avg[has])
            dmin[has] = np.minimum(dmin[has], dd)
        better = dmin < best
        best[better] = dmin[better]
        best_face[better] = f_local[better]
    d_plane = np.abs(np.einsum("ij,ij->i", cen - fc[best_face], nhat[best_face]))
    return dict(d_fc=d_fc[:, 0], d_exact=best, d_plane=d_plane, nearest_wall_face=wf[best_face])


def main():
    case_dir, stem, cache = cff.args()
    m = cff.load(case_dir, stem, cache)
    d = wall_distances(m)
    np.savez(os.path.join(cache, stem + "_walldist.npz"), **d)
    ref = cff.cell_field(m["cells"], "SV_WALL_DIST", m["n_cells"])
    # Wall-adjacent cells: c0 of wall faces.
    wf = wall_faces(m)
    adj = np.zeros(m["n_cells"], bool)
    adj[m["c0"][wf]] = True
    print("cells", ref.size, "wall-adjacent", adj.sum())
    for key in ("d_fc", "d_exact", "d_plane"):
        r = d[key] / ref - 1
        for label, mask in (("all", slice(None)), ("wall-adjacent", adj), ("interior", ~adj)):
            a = np.abs(r[mask])
            print(f"{key:8s} {label:14s} |d/SV_WALL_DIST-1|: median {np.median(a):.2e} p95 {np.quantile(a, .95):.2e} max {a.max():.2e}  signed median {np.median(r[mask]):+.2e}")


if __name__ == "__main__":
    main()
