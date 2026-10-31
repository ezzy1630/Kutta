"""Documented SST + intermittency source terms in double precision, shared by the
independent verification solvers in this directory. Nothing here is taken from
Kutta's kernels; every line cites its source.

Sources:
  TG    Ansys Fluent 2026 R1 Theory Guide, 4.8.2 (gamma equation, Eqs. 4-215..217
        and the unnumbered P_gamma, E_gamma, Tu_L, F_PG, Re_theta_c and constants)
        and 4.7.2.2 (k coupling, Eqs. 4-205..207), referenced by 4.8.3.
  TMR   NASA Turbulence Modeling Resource, "SST-2003-Menter-gamma-2015" page
        (Menter, Smirnov, Liu, Avancha 2015, FTaC 95:583): P_k^lim, F_on^lim,
        Re_theta_c^lim = 1100, F1 = max(F1_SST, F3), F3, R_y, gamma BCs.
  SST   Menter 1994 constants with Fluent's CDkw floor 1e-10 and production
        limiter 10 beta* k omega (the same base for both runs, so only the
        intermittency terms differ between them).
All arrays are wall units or SI with rho = 1; nu is the kinematic viscosity.
"""
import numpy as np

BSTAR, KAPPA, A1 = 0.09, 0.41, 0.31
SK1, SW1, B1 = 0.85, 0.5, 0.075
SK2, SW2, B2 = 1.0, 0.856, 0.0828

# TG 4.8.2 defaults.
CTU1, CTU2, CTU3 = 100.0, 1000.0, 1.0
CPG1, CPG2, CPG3, CPG1_LIM, CPG2_LIM = 1.0, -0.5, 0.0, 1.5, 3.0
F_LENGTH, CA2, CE2, SIGMA_GAMMA = 100.0, 0.06, 50.0, 1.0
# TMR: P_k^lim constants.
CK, CSEP, RE_THETA_C_LIM = 1.0, 1.0, 1100.0


def sst_blend(k, w, y, nu, dkdy, dwdy, S, gamma_model):
    """F1 (with TMR's F1 = max(F1, F3) when the model is on), F2, nu_t and blended constants."""
    cd = np.maximum(2 * SW2 / w * dkdy * dwdy, 1e-10)
    arg1 = np.minimum(np.maximum(np.sqrt(k) / (BSTAR * w * y), 500 * nu / (y**2 * w)), 4 * SW2 * k / (cd * y**2))
    F1 = np.tanh(arg1**4)
    if gamma_model:
        F1 = np.maximum(F1, np.exp(-((y * np.sqrt(k) / nu / 120) ** 8)))  # TMR F3, R_y = d sqrt(k)/nu
    arg2 = np.maximum(2 * np.sqrt(k) / (BSTAR * w * y), 500 * nu / (y**2 * w))
    F2 = np.tanh(arg2**2)
    nut = A1 * k / np.maximum(A1 * w, S * F2)
    blend = lambda a, b: F1 * a + (1 - F1) * b  # noqa: E731
    sk, sw, beta = blend(SK1, SK2), blend(SW1, SW2), blend(B1, B2)
    alpha = beta / BSTAR - sw * KAPPA**2 / np.sqrt(BSTAR)
    return dict(F1=F1, F2=F2, nut=nut, sk=sk, sw=sw, beta=beta, alpha=alpha)


def onset(k, w, y, nu, S, dvdy):
    """TG 4.8.2: Tu_L, lambda_thetaL (4-216, bounded by 4-217), F_PG, Re_theta_c, F_onset (4-215), F_turb."""
    RT = k / (nu * w)
    TuL = np.minimum(100 * np.sqrt(2 * k / 3) / (w * y), 100)
    lam = np.clip(-0.1111 * dvdy * y**2 / nu + 0.1875, -10, 10)
    FPG = np.where(lam >= 0, np.minimum(1 + CPG1 * lam, CPG1_LIM),
                   np.minimum(1 + CPG2 * lam + CPG3 * np.minimum(lam + 1, 0), CPG2_LIM))
    FPG = np.maximum(FPG, 0)  # TMR robustness limiter; inactive with C_PG3 = 0
    Rec = CTU1 + CTU2 * np.exp(-CTU3 * TuL * FPG)
    ReV = y**2 * S / nu
    Fon = np.maximum(np.minimum(ReV / (2.2 * Rec), 2) - np.maximum(1 - (RT / 3.5) ** 3, 0), 0)
    Fturb = np.exp(-((RT / 2) ** 4))
    return dict(RT=RT, TuL=TuL, lam=lam, FPG=FPG, Rec=Rec, ReV=ReV, Fon=Fon, Fturb=Fturb)


def gamma_rates(S, Omega, o):
    """P_gamma - E_gamma = a g (1 - g) - b g (c_e2 g - 1), TG 4.8.2."""
    return F_LENGTH * S * o["Fon"], CA2 * Omega * o["Fturb"]


def k_coupling(g, Pk, nut, S, Omega, y, nu, floor=0.1):
    """TG 4-206/207 (G*_k = gamma G~_k, Y*_k = min(max(gamma, 0.1), 1) Y_k) plus TMR's P_k^lim.
    Returns (production, destruction multiplier). `floor` exists only for ablation."""
    ReV = y**2 * S / nu
    Fonlim = np.minimum(np.maximum(ReV / (2.2 * RE_THETA_C_LIM) - 1, 0), 3)
    Plim = 5 * CK * np.maximum(g - 0.2, 0) * (1 - g) * Fonlim * np.maximum(3 * CSEP * nu - nut, 0) * S * Omega
    return g * Pk + Plim, np.minimum(np.maximum(g, floor), 1.0)
