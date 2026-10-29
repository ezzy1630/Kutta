"""Decode Fluent's k-omega near-wall treatment from saved wall fields.

Tests, per wall face / wall-adjacent cell:
  * what SV_WALL_SHEAR is (stress or force) and the tau_w law,
  * how SV_WALL_YPLUS and SV_WALL_YPLUS_UTAU are defined (distance, velocity scale),
  * the omega value in wall-adjacent cells and on wall faces,
  * the wall k boundary condition (face k vs cell k).
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.optimize import least_squares

import cff
import walldist

KAPPA, E_LOG, CMU, BSTAR = 0.4187, 9.793, 0.09, 0.09


def q(x):
    a = np.abs(x - 1)
    return f"ratio median {np.median(x):.6f} | |ratio-1| median {np.median(a):.2e} p95 {np.quantile(a, .95):.2e} max {a.max():.2e}"


def omega_wall(rho, mu, ustar, yplus, beta_i=0.072, c_calib=1 / 3, c_exp=1.3, kappa_sqrt_bstar=KAPPA * np.sqrt(BSTAR)):
    """TG Eq. 4-423..4-426 (log-law du+/dy+ = 1/(kappa y+))."""
    w_lam = 6 * c_calib / (beta_i * yplus**2)
    w_turb = 1 / (kappa_sqrt_bstar * yplus)
    return rho * ustar**2 / mu * w_lam * (1 + (w_turb / w_lam) ** c_exp) ** (1 / c_exp)


def main():
    case_dir, stem, cache = cff.args()
    m = cff.load(case_dir, stem, cache)
    n = m["n_cells"]
    C, F = m["cells"], m["faces"]
    nbf = max(z["hi"] for z in m["zones"] if z["type"] != 2)
    wf = walldist.wall_faces(m)
    c = m["c0"][wf]
    g = lambda name: cff.cell_field(C, name, n)[c]
    fw = lambda name: cff.cell_field(F, name, nbf)[wf]
    k, w, rho, mu, ydist = g("SV_K"), g("SV_O"), g("SV_DENSITY"), g("SV_MU_LAM"), g("SV_WALL_DIST")
    U = np.stack([g("SV_U"), g("SV_V"), g("SV_W")], 1)
    A = m["farea"][wf]
    area = np.linalg.norm(A, axis=1)
    nh = A / area[:, None]  # RH normal points into c0, i.e. into the fluid
    dvec = m["cen"][c] - m["fcent"][wf]
    d_n = np.einsum("ij,ij->i", dvec, nh)  # wall-normal centroid distance
    d_c = np.linalg.norm(dvec, axis=1)  # centroid-to-face-centroid distance
    Ut = U - np.einsum("ij,ij->i", U, nh)[:, None] * nh  # wall is stationary here
    ut = np.linalg.norm(Ut, axis=1)
    shear = fw("SV_WALL_SHEAR")
    yp_saved, ypu_saved = fw("SV_WALL_YPLUS"), fw("SV_WALL_YPLUS_UTAU")
    k_face, w_face = fw("SV_K"), fw("SV_O")
    assert np.all(fw("SV_U") == 0) and np.all(fw("SV_V") == 0) and np.all(fw("SV_W") == 0)

    print(f"wall faces {wf.size}; d_n>0 everywhere: {np.all(d_n > 0)}; d_n/d_c median {np.median(d_n / d_c):.5f} min {np.min(d_n / d_c):.3f}")
    print("\n1) SV_WALL_SHEAR")
    print(f"   |n.shear|/|shear| max {np.max(np.abs(np.einsum('ij,ij->i', shear, nh)) / np.linalg.norm(shear, axis=1)):.1e} (tangential)")
    print(f"   cos(shear, U_t) median {np.median(np.einsum('ij,ij->i', shear, Ut) / np.linalg.norm(shear, axis=1) / ut):+.6f}")
    print("   |shear| / (mu U_t/d_n)            :", q(np.linalg.norm(shear, axis=1) / (mu * ut / d_n)))
    tau = np.linalg.norm(shear, axis=1) / area
    print("   (|shear|/area) / (mu U_t/d_n)     :", q(tau / (mu * ut / d_n)))

    print("\n2) y+ definitions (candidates: distance d_n | d_c | SV_WALL_DIST; velocity scale u_tau | C_mu^.25 k^.5 | u* of Eq. 4-418)")
    u_tau = np.sqrt(tau / rho)
    u_k = CMU**0.25 * np.sqrt(k)
    for dname, d in (("d_n", d_n), ("d_c", d_c), ("SV_WALL_DIST", ydist)):
        u_star = np.sqrt((mu * ut / d + np.sqrt(CMU) * rho * k) / rho)
        for uname, u in (("u_tau", u_tau), ("C_mu^.25 k^.5", u_k), ("u*(4-418)", u_star)):
            print(f"   {dname:12s} {uname:14s} SV_WALL_YPLUS: {q(yp_saved / (rho * u * d / mu))}")
            print(f"   {'':12s} {'':14s} SV_WALL_YPLUS_UTAU: {q(ypu_saved / (rho * u * d / mu))}")

    d = d_n
    u_star = np.sqrt((mu * ut / d + np.sqrt(CMU) * rho * k) / rho)
    yplus = rho * u_star * d / mu
    print(f"\n   y+ range in this case: {yplus.min():.2f} .. {yplus.max():.2f} (median {np.median(yplus):.2f})")

    print("\n3) tau_w law, Eq. 4-413..4-417 with d_n and u* (4-418)")

    def tau_model(kappa, E):
        up_lam = yplus
        up_turb = np.log(E * np.maximum(yplus, 0.2)) / kappa
        up = (up_lam**-4 + up_turb**-4) ** -0.25
        return rho * u_star * ut / up

    print("   documented kappa=0.4187, E=9.793:", q(tau / tau_model(KAPPA, E_LOG)))
    fit = least_squares(lambda p: np.log(tau / tau_model(*p)), [0.41, 9.0])
    print(f"   free fit kappa={fit.x[0]:.4f} E={fit.x[1]:.3f} rms log err {np.sqrt(np.mean(fit.fun**2)):.2e}")

    print("\n4) omega in wall-adjacent cells (SV_O cell) vs candidates")
    nu = mu / rho
    cands = {
        "Eq 4-423, beta_i=0.072 (std k-w), C_calib=1/3, C_exp=1.3": omega_wall(rho, mu, u_star, yplus, 0.072),
        "Eq 4-423, beta_i=0.075 (SST beta_i1)": omega_wall(rho, mu, u_star, yplus, 0.075),
        "Eq 4-423, beta_i=0.0828 (SST beta_i2)": omega_wall(rho, mu, u_star, yplus, 0.0828),
        "viscous 6 nu/(0.075 y^2)": 6 * nu / (0.075 * d**2),
        "TMR wall 60 nu/(0.075 y^2)": 60 * nu / (0.075 * d**2),
        "log u*/(sqrt(beta*) kappa y)": u_star / (np.sqrt(BSTAR) * KAPPA * d),
        "Menter automatic sqrt(visc^2 + log^2)": np.sqrt((6 * nu / (0.075 * d**2)) ** 2 + (u_star / (np.sqrt(BSTAR) * KAPPA * d)) ** 2),
    }
    for name, val in cands.items():
        print(f"   {name:58s} {q(w / val)}")
    fit = least_squares(
        lambda p: np.log(w / omega_wall(rho, mu, u_star, yplus, beta_i=p[0], c_exp=p[1], kappa_sqrt_bstar=p[2])),
        [0.075, 1.3, 0.12],
    )
    print(f"   free fit (C_calib=1/3 fixed): beta_i={fit.x[0]:.6f} C_exp={fit.x[1]:.5f} kappa*sqrt(beta*)={fit.x[2]:.6f} "
          f"(0.4187*0.3={KAPPA * 0.3:.6f})  rms log err {np.sqrt(np.mean(fit.fun**2)):.2e} max {np.abs(fit.fun).max():.2e}")
    print("   same with d = SV_WALL_DIST instead of d_n:")
    u2 = np.sqrt((mu * ut / ydist + np.sqrt(CMU) * rho * k) / rho)
    print(f"   {'':58s} {q(w / omega_wall(rho, mu, u2, rho * u2 * ydist / mu, 0.072))}")
    print("   same with u* = C_mu^.25 k^.5 only:")
    print(f"   {'':58s} {q(w / omega_wall(rho, mu, u_k, rho * u_k * d / mu, 0.072))}")

    print("\n5) wall-face SV_O")
    print("   face/cell:", q(w_face / w))
    fit = least_squares(
        lambda p: np.log(w_face / (rho * u_star**2 / mu * (p[0] / yplus**2) * (1 + ((p[1] / yplus) / (p[0] / yplus**2)) ** p[2]) ** (1 / p[2]))),
        [27.8, 8.0, 1.3],
    )
    print(f"   free fit w+ = blend(a/y+^2, b/y+, c): a={fit.x[0]:.3f} (6/0.072={6 / 0.072:.3f}) b={fit.x[1]:.3f} "
          f"(2/(0.4187*0.3)={2 / (KAPPA * 0.3):.3f}) c={fit.x[2]:.4f}  rms log err {np.sqrt(np.mean(fit.fun**2)):.2e}")
    fixed = omega_wall(rho, mu, u_star, yplus, beta_i=0.072, c_calib=1.0, kappa_sqrt_bstar=KAPPA * 0.3 / 2)
    print("   fixed a=6/0.072, b=2/(kappa sqrt(beta*)), c=1.3:", q(w_face / fixed))

    print("\n6) wall k: face SV_K / cell SV_K:", q(k_face / k))
    mut = g("SV_MU_T")
    print("   wall-cell mu_t / (rho k/omega):", q(mut / (rho * k / w)))


if __name__ == "__main__":
    main()
