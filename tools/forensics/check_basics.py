"""Sanity checks: geometry closure, orientation, field ranges, convergence."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import h5py
import numpy as np
import cff

case_dir, stem, cache = cff.args()
m = cff.load(case_dir, stem, cache)
n = m["n_cells"]; c0, c1 = m["c0"], m["c1"]; A = m["farea"]
inn = c1 >= 0
clos = np.zeros((n, 3)); np.add.at(clos, c0, -A); np.add.at(clos, c1[inn], A[inn])
print("cells", n, "faces", c0.size, "vol min/sum", m["vol"].min(), m["vol"].sum())
print("max |closure|/area^(..)", np.max(np.linalg.norm(clos, axis=1) / m["vol"] ** (2 / 3)))
# orientation: RH normal points into c0  <=> (cen[c0]-fcent).A > 0
d = np.einsum("ij,ij->i", m["cen"][c0] - m["fcent"], A)
print("faces with normal into c0:", np.mean(d > 0))
for z in m["zones"]: print(z)
C = m["cells"]
for k in sorted(C):
    v = C[k][1]; print(k, v.shape, np.nanmin(v, axis=0), np.nanmax(v, axis=0))
F = m["faces"]
for k in sorted(F):
    ids, v = F[k]; print("face", k, ids.min(), ids.max(), v.shape)
flux = cff.cell_field(F, "SV_FLUX", c0.size)
for z in m["zones"]: print(z["name"], flux[z["lo"]:z["hi"]].sum())
res = np.zeros(n); np.add.at(res, c0, flux); np.add.at(res, c1[inn], -flux[inn])
print("cell mass imbalance / max |flux| :", np.abs(res).max() / np.abs(flux).max())
with h5py.File(os.path.join(case_dir, stem + ".dat.h5")) as f:
    r = f["results/residuals"]
    def v(name, o):
        if isinstance(o, h5py.Dataset) and name.endswith("data"):
            d = o[:]; it = r[name.replace("data", "iterations")][:]
            print(name, d.shape, "iters", it[0], it[-1], "last raw/scale", d[-1, 0] / d[-1, 1] if d[-1, 1] else d[-1])
    r.visititems(v)
    print(f["settings/Data Variables"][0][:600])
