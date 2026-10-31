"""Zero-pressure-gradient flat plate, SST and SST + intermittency, from the documented equations alone.

Parabolic boundary-layer equations in double precision: U_e = rho = 1, lengths in
metres, nu = 1/Re_unit. Implicit upwind marching in x, conservative central
diffusion in y, Picard iterations per station. Boundary-layer approximation:
S = Omega = |du/dy| and dV/dy = dv/dy = -du/dx (TMR's n . grad(n . U) on a flat
plate). Wall: u = 0, k = 0, omega = 60 nu/(beta1 dy1^2) (Menter's wall value);
gamma zero flux (TG 4.7.2.1 text, TMR). Top: freestream SST decay with F1 = 0,
U dk/dx = -beta* k w, U dw/dx = -beta2 w^2, gamma = 1. The freestream starts
x_in upstream of the leading edge with the inlet Tu and mu_t/mu that Kutta uses.

  build-py/bin/python tools/transition/verify/boundary_layer.py t3a [--out DIR] [--check]

writes DIR/independent_<case>.csv (x, Re_x, Cf SST, Cf gamma). --check reruns the
case on finer grids and smaller steps and prints the changes.
"""
import argparse
import pathlib

import numpy as np
from scipy.linalg import solve_banded

import model as m

# Re_unit, inlet Tu, inlet mu_t/mu, inlet distance upstream of the leading edge, plate end.
# T3 values are tools/transition/flatplate.py's fits (validation/results/transition/runs/*_fit.json);
# TMR is bench_tmr's high-tu mode on the 2DZP plate (inflow at x = -0.33333).
CASES = {
    "t3a": (343807.450411224, 0.046778009427568656, 11.926800016097618, 0.1, 1.7),
    "t3b": (640133.779264214, 0.07420386830418466, 100.23209424690732, 0.1, 1.7),
    "t3am": (1305968.2639930528, 0.010211100138978212, 8.057295339838335, 0.1, 1.7),
    "tmr": (5e6, 0.05, 100.0, 0.33333, 2.0),
}


def stretched(n, first, length):
    lo, hi = 1.0 + 1e-12, 2.0
    for _ in range(200):
        r = 0.5 * (lo + hi)
        lo, hi = (r, hi) if first * (r**n - 1) / (r - 1) < length else (lo, r)
    d = first * r ** np.arange(n)
    return np.concatenate([[0.0], np.cumsum(d * length / d.sum())])


def run(re_unit, tu_in, ratio, x_in, x_end, gamma_model, H=0.12, ny=260, dy1=None,
        dx0=2e-6, dxmax=2e-3, growth=1.01, picard=4, profile_x=(), edge=None):
    """Returns x, Cf and {x: profile} for the requested stations. `edge` is an
    optional edge velocity U_e(x) (default 1, zero pressure gradient); the
    pressure gradient then enters as U_e dU_e/dx."""
    nu = 1.0 / re_unit
    y = stretched(ny, dy1 or 0.2 / (re_unit * np.sqrt(0.003)), H)
    n = len(y)
    k0, dyf = 1.5 * tu_in**2, np.diff(y)
    w0 = k0 / (nu * ratio)
    dyc = np.zeros(n); dyc[1:-1] = 0.5 * (y[2:] - y[:-2])
    yy = np.maximum(y, 1e-30)
    w_wall = 60 * nu / (m.B1 * y[1] ** 2)

    def freestream(x):
        t = x + x_in
        return k0 * (1 + m.B2 * w0 * t) ** (-m.BSTAR / m.B2), w0 / (1 + m.B2 * w0 * t)

    def ddy(q):
        d = np.empty(n)
        d[1:-1] = (q[2:] - q[:-2]) / (y[2:] - y[:-2])
        d[0], d[-1] = (q[1] - q[0]) / dyf[0], (q[-1] - q[-2]) / dyf[-1]
        return d

    def march(old, G, src, sink, u, v, dx, wall, top, zero_flux=False):
        """u (phi - old)/dx + v dphi/dy = d/dy(G dphi/dy) + src - sink phi, upwind in v."""
        Gf = 0.5 * (G[1:] + G[:-1])
        i = np.arange(1, n - 1)
        aW, aE = np.zeros(n), np.zeros(n)
        aW[i] = Gf[i - 1] / dyf[i - 1] / dyc[i] + np.maximum(v[i], 0) / dyf[i - 1]
        aE[i] = Gf[i] / dyf[i] / dyc[i] - np.minimum(v[i], 0) / dyf[i]
        diag, b = np.ones(n), np.zeros(n)
        diag[i] = aW[i] + aE[i] + u[i] / dx + sink[i]
        b[i] = u[i] / dx * old[i] + src[i]
        lower, upper = -aW, -aE
        upper[0], b[0] = (-1.0, 0.0) if zero_flux else (0.0, wall)
        lower[-1], b[-1] = 0.0, top
        ab = np.zeros((3, n)); ab[0, 1:] = upper[:-1]; ab[1] = diag; ab[2, :-1] = lower[1:]
        return solve_banded((1, 1), ab, b)

    x = 1e-7
    ke, we = freestream(x)
    u, v = np.ones(n), np.zeros(n)
    u[0] = 0.0
    k, w, g = np.full(n, ke), np.full(n, we), np.ones(n)
    k[0], w[0] = 0.0, w_wall
    xs, cfs, profiles, pending = [], [], {}, sorted(profile_x)
    dx = dx0
    while x < x_end:
        xn = min(x + dx, x_end)
        dx = xn - x
        ke, we = freestream(xn)
        ue = 1.0 if edge is None else edge(xn)
        dpdx = 0.0 if edge is None else ue * (ue - edge(x)) / dx  # -dp/dx = U_e dU_e/dx
        uo, ko, wo, go = u.copy(), k.copy(), w.copy(), g.copy()
        for _ in range(picard):
            S = np.abs(ddy(u))
            kk = np.maximum(k, 1e-20)
            b = m.sst_blend(kk, w, yy, nu, ddy(k), ddy(w), S, gamma_model)
            b["F1"][0] = 1.0
            nut = b["nut"]; nut[0] = 0.0
            u = march(uo, nu + nut, np.full(n, dpdx), np.zeros(n), u, v, dx, 0.0, ue)
            dudx = (u - uo) / dx
            v = np.concatenate([[0.0], -np.cumsum(0.5 * (dudx[1:] + dudx[:-1]) * dyf)])
            S = np.abs(ddy(u))
            Pk = np.minimum(nut * S * S, 10 * m.BSTAR * kk * w)
            prod, dmul = Pk, 1.0
            if gamma_model:
                prod, dmul = m.k_coupling(g, Pk, nut, S, S, yy, nu)
            dk, dw = ddy(k), ddy(w)
            k = np.maximum(march(ko, nu + b["sk"] * nut, prod, dmul * m.BSTAR * w, u, v, dx, 0.0, ke), 1e-20)
            cross = 2 * (1 - b["F1"]) * m.SW2 / w * dk * dw
            src = b["alpha"] * S * S + np.maximum(cross, 0) + b["beta"] * w * w
            sink = 2 * b["beta"] * w + np.maximum(-cross, 0) / w
            w = np.maximum(march(wo, nu + b["sw"] * nut, src, sink, u, v, dx, w_wall, we), 1e-20)
            if gamma_model:
                kk = np.maximum(k, 1e-20)
                nut = m.A1 * kk / np.maximum(m.A1 * w, S * b["F2"]); nut[0] = 0.0
                o = m.onset(kk, w, yy, nu, S, -dudx)
                a, bb = m.gamma_rates(S, S, o)
                # Same steady source as P - E; positive part explicit, quadratic part implicit.
                g = np.clip(march(go, nu + nut / m.SIGMA_GAMMA, (a + bb) * g, (a + m.CE2 * bb) * g, u, v, dx,
                                  0.0, 1.0, zero_flux=True), 1e-7, 1.0)
        x = xn
        xs.append(x); cfs.append(2 * nu * u[1] / y[1] / u[-1] ** 2)  # on the local U_e, as Kutta's T3 Cf
        while pending and x >= pending[0]:
            profiles[pending.pop(0)] = dict(y=y.copy(), u=u.copy(), k=k.copy(), w=w.copy(), g=g.copy(),
                                            nut=nut.copy(), cf=cfs[-1], nu=nu)
        dx = min(dx * growth, dxmax)
    return np.array(xs), np.array(cfs), profiles


def transition(rex, cf):
    """Onset (Cf minimum before the rise) and end (first maximum after it) of a Cf(Re_x) curve."""
    keep = rex > 3e3
    r, c = rex[keep], cf[keep]
    rise = np.flatnonzero(c > 1.1 * np.minimum.accumulate(c))
    if len(rise) == 0:
        return np.nan, np.nan
    i = int(np.argmin(c[:rise[0] + 1]))
    j = rise[0]
    while j + 1 < len(c) and c[j + 1] >= c[j]:
        j += 1
    return r[i], r[j]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case", choices=list(CASES))
    ap.add_argument("--out", default=str(pathlib.Path(__file__).resolve().parents[3] / "validation/results/transition_verification"))
    ap.add_argument("--check", action="store_true", help="print the grid/step sensitivity instead of writing")
    args = ap.parse_args()
    c = CASES[args.case]
    if args.check:
        base = None
        for name, kw in [("base", {}), ("ny 400", dict(ny=400)), ("dy1 / 2", dict(dy1=0.1 / (c[0] * np.sqrt(0.003)))),
                         ("dx / 4", dict(dxmax=5e-4, dx0=5e-7)), ("picard 8", dict(picard=8)), ("H 0.2", dict(H=0.2, ny=300))]:
            xg, cg, _ = run(*c, True, **kw)
            xs, cs, _ = run(*c, False, **kw)
            on, end = transition(xg * c[0], cg)
            probe = np.geomspace(2e-3, c[4], 60)
            if base is None:  # regions fixed by the base run
                lam, turb = probe * c[0] < 0.8 * on, probe * c[0] > 1.5 * end
            vals = dict(gl=np.interp(probe[lam], xg, cg), gt=np.interp(probe[turb], xg, cg), s=np.interp(probe, xs, cs))
            if base is None:
                base, base_on = vals, on
            dev = {key: 100 * np.max(abs(vals[key] / base[key] - 1), initial=0) for key in vals}
            print(f"{name:9s} onset Re_x {on:.4e} ({100*(on/base_on-1):+.2f}%), gamma laminar {dev['gl']:.2f}%, "
                  f"gamma turbulent {dev['gt']:.2f}%, SST {dev['s']:.2f}%")
        return
    xg, cg, _ = run(*c, True)
    xs, cs, _ = run(*c, False)
    out = pathlib.Path(args.out) / f"independent_{args.case}.csv"
    with open(out, "w") as f:
        f.write("x,re_x,cf_sst,cf_gamma\n")
        for x in np.unique(np.r_[np.geomspace(1e-3, c[4], 400), c[4]]):
            f.write(f"{x:.6e},{x*c[0]:.6e},{np.interp(x, xs, cs):.6e},{np.interp(x, xg, cg):.6e}\n")
    print(out)


if __name__ == "__main__":
    main()
