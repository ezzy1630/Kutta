"""Minimal Fluent CFF (.cas.h5/.dat.h5) reader plus finite-volume geometry.

Parse-only: reads arrays with h5py, never evaluates anything from the file.
All ids in the file are 1-based; everything returned here is 0-based.

Forensics scripts in this directory (see docs/research/fluent-sst-forensics.md):
    python -I <script>.py [case_dir] [stem] [cache_dir]
Defaults: validation/fluent_examples/mixing_elbow, mixing_elbow, /tmp/kutta-forensics/cache.
Order: settings.py, check_basics.py, walldist.py, mut_limiter.py, wall_treatment.py, residuals.py;
then strain_gradient.py, face_flux.py (fluent-flux-gradient-forensics.md) and momentum_residual.py,
energy_residual.py, turbulence_residual.py (fluent-momentum-forensics.md; building blocks in momentum.py).
Needs numpy, scipy, h5py (e.g. a uv venv outside the repo).
"""
import os
import sys

import h5py
import numpy as np


def _sections(group):
    """Yield (minId, maxId, dataset_or_group) for each numbered section."""
    for key in sorted(group.keys(), key=lambda k: int(k) if k.isdigit() else -1):
        if not key.isdigit():
            continue
        obj = group[key]
        yield int(obj.attrs["minId"][0]), int(obj.attrs["maxId"][0]), obj


def read_mesh(cas_path):
    with h5py.File(cas_path, "r") as f:
        m = f["meshes/1"]
        n_nodes = int(m.attrs["nodeCount"][0])
        n_cells = int(m.attrs["cellCount"][0])
        n_faces = int(m.attrs["faceCount"][0])
        nodes = np.zeros((n_nodes, 3))
        for lo, hi, ds in _sections(m["nodes/coords"]):
            nodes[lo - 1 : hi] = ds[:]
        c0 = np.zeros(n_faces, dtype=np.int64)
        c1 = np.full(n_faces, -1, dtype=np.int64)
        for lo, hi, ds in _sections(m["faces/c0"]):
            c0[lo - 1 : hi] = ds[:].astype(np.int64) - 1
        for lo, hi, ds in _sections(m["faces/c1"]):
            c1[lo - 1 : hi] = ds[:].astype(np.int64) - 1
        nnodes = np.zeros(n_faces, dtype=np.int64)
        flat = []
        for lo, hi, grp in _sections(m["faces/nodes"]):
            nnodes[lo - 1 : hi] = grp["nnodes"][:]
            flat.append(grp["nodes"][:].astype(np.int64) - 1)
        flat = np.concatenate(flat)
        assert flat.size == nnodes.sum()
        zt = m["faces/zoneTopology"]
        names = zt["name"][0].decode().split(";")
        zones = []
        for i, zid in enumerate(zt["id"][:]):
            zones.append(
                dict(
                    id=int(zid),
                    name=names[i],
                    type=int(zt["zoneType"][i]),
                    lo=int(zt["minId"][i]) - 1,
                    hi=int(zt["maxId"][i]),
                )
            )
    return dict(nodes=nodes, c0=c0, c1=c1, nnodes=nnodes, fnodes=flat, zones=zones, n_cells=n_cells)


def face_geometry(mesh):
    """Face area vector (RH rule from node order) and centroid.

    Polygon split into triangles about the node average; centroid is the
    area-weighted triangle centroid (standard FV practice, also what OpenFOAM does).
    """
    nodes, nn, fn = mesh["nodes"], mesh["nnodes"], mesh["fnodes"]
    nf = nn.size
    ptr = np.concatenate([[0], np.cumsum(nn)])
    face_of = np.repeat(np.arange(nf), nn)
    pts = nodes[fn]
    avg = np.zeros((nf, 3))
    np.add.at(avg, face_of, pts)
    avg /= nn[:, None]
    nxt = np.arange(fn.size) + 1
    nxt[ptr[1:] - 1] = ptr[:-1]
    a = pts
    b = nodes[fn[nxt]]
    c = avg[face_of]
    tri_area = 0.5 * np.cross(a - c, b - c)
    tri_cent = (a + b + c) / 3.0
    area = np.zeros((nf, 3))
    np.add.at(area, face_of, tri_area)
    # Weight triangle centroids by projection of their area on the face normal.
    nhat = area / np.linalg.norm(area, axis=1)[:, None]
    w = np.einsum("ij,ij->i", tri_area, nhat[face_of])
    cent = np.zeros((nf, 3))
    np.add.at(cent, face_of, tri_cent * w[:, None])
    wsum = np.zeros(nf)
    np.add.at(wsum, face_of, w)
    cent /= wsum[:, None]
    return area, cent


def cell_geometry(mesh, farea, fcent):
    """Cell volume and centroid by pyramid decomposition about the face-centroid average."""
    c0, c1, n = mesh["c0"], mesh["c1"], mesh["n_cells"]
    interior = c1 >= 0
    cnt = np.bincount(c0, minlength=n) + np.bincount(c1[interior], minlength=n)
    est = np.zeros((n, 3))
    np.add.at(est, c0, fcent)
    np.add.at(est, c1[interior], fcent[interior])
    est /= cnt[:, None]
    # RH normal points INTO c0, so outward area for c0 is -farea, for c1 is +farea.
    vol = np.zeros(n)
    cen = np.zeros((n, 3))
    for cells, sign, mask in ((c0, -1.0, slice(None)), (c1, 1.0, interior)):
        cc = cells[mask]
        h = np.einsum("ij,ij->i", fcent[mask] - est[cc], sign * farea[mask])
        pv = h / 3.0
        pc = 0.75 * fcent[mask] + 0.25 * est[cc]
        np.add.at(vol, cc, pv)
        np.add.at(cen, cc, pc * pv[:, None])
    cen /= vol[:, None]
    return vol, cen, cnt


def read_fields(dat_path, where):
    """Return {name: (ids0, array)} for /results/1/phase-1/<where>."""
    out = {}
    with h5py.File(dat_path, "r") as f:
        g = f["results/1/phase-1/" + where]
        for name in g.keys():
            if not name.startswith("SV_"):
                continue
            ids, vals = [], []
            for lo, hi, ds in _sections(g[name]):
                ids.append(np.arange(lo - 1, hi))
                vals.append(ds[:])
            out[name] = (np.concatenate(ids), np.concatenate(vals))
    return out


def cell_field(fields, name, n):
    ids, v = fields[name]
    arr = np.full((n,) + v.shape[1:], np.nan)
    arr[ids] = v
    return arr


def load(case_dir, stem, cache_dir):
    """Load mesh, geometry and fields; cache derived geometry in cache_dir."""
    cas = os.path.join(case_dir, stem + ".cas.h5")
    dat = os.path.join(case_dir, stem + ".dat.h5")
    mesh = read_mesh(cas)
    cache = os.path.join(cache_dir, stem + "_geom.npz")
    if os.path.exists(cache):
        z = np.load(cache)
        farea, fcent, vol, cen, cnt = (z[k] for k in ("farea", "fcent", "vol", "cen", "cnt"))
    else:
        farea, fcent = face_geometry(mesh)
        vol, cen, cnt = cell_geometry(mesh, farea, fcent)
        np.savez(cache, farea=farea, fcent=fcent, vol=vol, cen=cen, cnt=cnt)
    mesh.update(farea=farea, fcent=fcent, vol=vol, cen=cen, nfaces_per_cell=cnt)
    mesh["cells"] = read_fields(dat, "cells")
    mesh["faces"] = read_fields(dat, "faces")
    return mesh


def args():
    repo = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    case_dir = sys.argv[1] if len(sys.argv) > 1 else os.path.join(repo, "validation/fluent_examples/mixing_elbow")
    stem = sys.argv[2] if len(sys.argv) > 2 else "mixing_elbow"
    cache = sys.argv[3] if len(sys.argv) > 3 else "/tmp/kutta-forensics/cache"
    os.makedirs(cache, exist_ok=True)
    return case_dir, stem, cache
