# /// script
# requires-python = "==3.14.*"
# dependencies = ["numpy==2.5.4"]
# ///
"""Small analytic wheel above a road for starter/CLI coverage, not aero accuracy.

uv run tools/fsae_aero/wheel_starter_mesh.py --out tests/data/meshes/fsae_wheel.msh
Centre (0.4, 0.25), radius 0.2 m; +x flow, +y up, +z rolling axis.
"""
import argparse
import math
from pathlib import Path

import numpy as np

import fluent_msh


def build(out):
    cx, cy, radius = 0.4, 0.25, 0.2
    bounds = (-1, 2, 0, 1.2)
    angles = sorted(set([2 * math.pi * i / 32 for i in range(32)] + [
        math.atan2(y - cy, x - cx) % (2 * math.pi)
        for x in bounds[:2] for y in bounds[2:]
    ]))
    n, rings = len(angles), 6
    points = []
    for j in range(rings + 1):
        for theta in angles:
            dx, dy = math.cos(theta), math.sin(theta)
            limits = [(x - cx) / dx for x in bounds[:2] if abs(dx) > 1e-14]
            limits += [(y - cy) / dy for y in bounds[2:] if abs(dy) > 1e-14]
            outer = min(t for t in limits if t > 0)
            r = radius + (outer - radius) * (j / rings) ** 1.5
            points.append([cx + r * dx, cy + r * dy, 0])
    points = np.array(points)
    cells = [[j*n+i, (j+1)*n+i, (j+1)*n+(i+1)%n, j*n+(i+1)%n]
             for j in range(rings) for i in range(n)]
    zones = {"inlet": [], "outlet": [], "ground": [], "top": [], "wheel_front": []}
    for i in range(n):
        zones["wheel_front"].append([i, (i+1)%n])
        edge = [rings*n+i, rings*n+(i+1)%n]
        x, y, _ = points[edge].mean(axis=0)
        name = min([("inlet", abs(x-bounds[0])), ("outlet", abs(x-bounds[1])),
                    ("ground", abs(y-bounds[2])), ("top", abs(y-bounds[3]))], key=lambda p: p[1])[0]
        zones[name].append(edge)
    return fluent_msh.write(out, points, [(3, np.array(cells))],
                            [(name, "wall", np.array(edges)) for name, edges in zones.items()],
                            "analytic 2D wheel starter: centre .4 .25, radius .2; no validation data")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    print(build(args.out))
