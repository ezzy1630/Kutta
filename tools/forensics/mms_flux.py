"""Isolate interpolation curvature in the V2 extruded-hex exact fields.

Run with the forensics venv: python -I tools/forensics/mms_flux.py.
This double-precision, equal-diagonal stencil probe is NOT a solver run: no
pressure correction, momentum solve, or limiter is used. It shows the local
boundary-divergence order, including the point formula's own boundary defect.
Actual converged solution/error layers are printed by test_navier_stokes.cpp.
"""
import numpy as np


def f(x):
    return x * x * (1 - x) ** 2


def df(x):
    return 2 * x - 6 * x * x + 4 * x * x * x


def wall_lsq(u, h, axis):
    """Exact two-row 1/d^2 LSQ, including the zero wall value at distance h/2."""
    v = np.moveaxis(u, axis, 0)
    g = np.empty_like(v)
    g[1:-1] = (v[2:] - v[:-2]) / (2 * h)
    g[0] = (v[1] + v[0]) / (2 * h)
    g[-1] = -(v[-2] + v[-1]) / (2 * h)
    return np.moveaxis(g, 0, axis)


def main():
    for n in (8, 16, 32, 64, 128):
        h = 1 / n
        x, y = np.meshgrid((np.arange(n) + .5) * h, (np.arange(n) + .5) * h, indexing="ij")
        u = np.stack((100 * f(x) * df(y), -100 * df(x) * f(y)), axis=-1)
        grad = [wall_lsq(u, h, axis) for axis in (0, 1)]
        fluxes = {name: [] for name in ("point", "side", "peclet")}
        for axis in (0, 1):
            v = np.moveaxis(u[..., axis], axis, 0)
            g = np.moveaxis(grad[axis][..., axis], axis, 0)
            point = .5 * (v[:-1] + v[1:])
            side = point + .25 * h * (g[:-1] - g[1:])
            # F and D divided by face area, with rho=1 and mu=.1 as in V2.
            blend = np.abs(side) / (np.abs(side) + 2 * .1 / h)
            variants = (("point", point), ("side", side), ("peclet", point + blend * (side - point)))
            for name, interior in variants:
                face = np.zeros((n + 1, n))
                face[1:-1] = interior
                fluxes[name].append(np.moveaxis(np.diff(face, axis=0) / h, 0, axis))
        layer = (np.minimum.reduce((x, 1 - x, y, 1 - y)) / h).astype(int)
        for name, components in fluxes.items():
            divergence = components[0] + components[1]
            print(n, name, "div_L1", np.mean(np.abs(divergence)),
                  "wall_abs_mean", np.mean(np.abs(divergence[layer == 0])),
                  "deep_abs_mean", np.mean(np.abs(divergence[layer >= 2])),
                  "wall_max", np.max(np.abs(divergence[layer == 0])))


if __name__ == "__main__":
    main()
