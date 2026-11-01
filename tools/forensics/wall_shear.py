"""Split Kutta's V6 wall-force errors face by face against Fluent's saved mixing elbow.

The V6 slow test gates each wall zone's integrated shear force. An integrated
force can agree while the faces do not, so this compares, per wall face, the
shear force on the wall, and per wall cell the tangential velocity and k, and
splits each zone into five bands along the coordinate it spans most. It also
reports k's volume-weighted bias and its error by wall layer.

Kutta's side comes from the V6 slow test with KUTTA_V6_WALL_DIR set, which
writes <backend>-<mode>-walls.csv (force on the wall per face, -shear*|A|,
and the wall cell in Fluent's order) and <backend>-<mode>-cells.csv.
Usage: python -I wall_shear.py <case_dir> <stem> <cache_dir> <wall_dir> <backend-mode>...
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.spatial import cKDTree

import cff
import fvops
import mut_limiter as ml


def ratio(a, b, w):
    return np.sum(a * w) / np.sum(b * w)


def main():
    case_dir, stem, cache = sys.argv[1:4]
    wall_dir, runs = sys.argv[4], sys.argv[5:]
    m = cff.load(case_dir, stem, cache)
    ctx = fvops.Ctx(m)
    g = lambda name: cff.cell_field(m["cells"], name, ctx.n)  # noqa: E731
    U = np.stack([g("SV_U"), g("SV_V"), g("SV_W")], 1)
    K = g("SV_K")
    # SV_WALL_SHEAR is the per-face force on the fluid; the wall receives its negative.
    F = -cff.cell_field(m["faces"], "SV_WALL_SHEAR", ctx.nbf)
    faces = cKDTree(ctx.fc[: ctx.nbf])
    layer = ml.wall_layers(m, 8)
    for run in runs:
        w = np.genfromtxt(os.path.join(wall_dir, run + "-walls.csv"), delimiter=",", names=True, dtype=None, encoding=None)
        c = np.genfromtxt(os.path.join(wall_dir, run + "-cells.csv"), delimiter=",", names=True)
        distance, face = faces.query(np.stack([w["x"], w["y"], w["z"]], 1))
        assert distance.max() < 1e-7 * np.abs(ctx.fc).max(), "wall faces do not match the case"
        T = np.stack([w["fx"], w["fy"], w["fz"]], 1)
        u = np.stack([c["u"], c["v"], c["w"]], 1)
        print(f"{run}")
        for patch in sorted(set(w["patch"])):
            s = w["patch"] == patch
            f, cell = face[s], w["cell"][s]
            n = ctx.nhat[f]
            tangential = lambda V: np.linalg.norm(V - np.einsum("ij,ij->i", V, n)[:, None] * n, axis=1)  # noqa: E731
            uk, uf = tangential(u[cell]), tangential(U[cell])
            tk, tf = np.linalg.norm(T[s], axis=1), np.linalg.norm(F[f], axis=1)
            e = np.linalg.norm(T[s] - F[f], axis=1)
            total = np.linalg.norm(T[s].sum(0) - F[f].sum(0)) / np.linalg.norm(F[f].sum(0))
            print(f"  {patch:12s} faces {s.sum():5d}: integrated force error {total * 100:.2f}% | per face L2 {np.sqrt((e**2).sum() / (tf**2).sum()) * 100:.2f}%"
                  f" L1 {e.sum() / tf.sum() * 100:.2f}% | wall-cell |U_t| L2 {np.sqrt(((uk - uf)**2).sum() / (uf**2).sum()) * 100:.2f}%"
                  f" | wall-cell k L2 {np.sqrt(((c['k'][cell] - K[cell])**2).sum() / (K[cell]**2).sum()) * 100:.1f}%")
            x = ctx.fc[f]
            axis = int(np.argmax(x.max(0) - x.min(0)))
            edges = np.quantile(x[:, axis], np.linspace(0, 1, 6))
            bands = []
            for lo, hi in zip(edges[:-1], edges[1:]):
                b = (x[:, axis] >= lo) & (x[:, axis] <= hi)
                bands.append(f"[{lo:+.3f},{hi:+.3f}] |tau| {ratio(tk[b], tf[b], 1):.3f} |U_t| {ratio(uk[b], uf[b], tf[b]):.3f} k {ratio(c['k'][cell][b], K[cell][b], tf[b]):.3f}")
            print(f"    Kutta/Fluent by {'xyz'[axis]} band: " + " | ".join(bands))
        total = (K**2).sum()
        by_layer = " ".join(f"L{L}:{np.sqrt(((c['k'] - K)[(layer == L) if L < 4 else (layer >= 4)]**2).sum() / total) * 100:.1f}" for L in (1, 2, 3, 4))
        print(f"  k: L2 {np.sqrt(((c['k'] - K)**2).sum() / total) * 100:.2f}% (by wall layer, L4 = deeper: {by_layer}),"
              f" volume-weighted bias {((c['k'] - K) * ctx.vol).sum() / (K * ctx.vol).sum() * 100:+.1f}%")


if __name__ == "__main__":
    main()
