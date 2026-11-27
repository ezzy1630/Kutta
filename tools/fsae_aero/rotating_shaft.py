"""Generate one grid in the workbench rotating-shaft verification family."""
import argparse
import math
from pathlib import Path

import numpy as np

import fluent_msh


def mesh(path, radial, angular):
    def node(ring, angle):
        return ring * angular + angle % angular

    points = []
    for j in range(radial + 1):
        radius = 0.1 + 0.1 * j / radial
        for i in range(angular):
            theta = 2 * math.pi * i / angular
            points.append([radius * math.cos(theta), radius * math.sin(theta), 0.])
    cells = [[node(j, i), node(j+1, i), node(j+1, i+1), node(j, i+1)]
             for j in range(radial) for i in range(angular)]
    walls = [(name, "wall", np.array([[node(j, i), node(j, i+1)] for i in range(angular)]))
             for name, j in [("shaft", 0), ("housing", radial)]]
    description = (f"Concentric circular walls: inner radius 0.1 m, outer radius 0.2 m; "
                   f"{radial} radial by {angular} angular cells.")
    return fluent_msh.write(path, np.array(points), [(3, np.array(cells))], walls, description)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--radial", type=int, default=32)
    parser.add_argument("--angular", type=int, default=144)
    args = parser.parse_args()
    if args.radial < 1 or args.angular < 3:
        parser.error("at least one radial cell and three circumferential cells are required")
    mesh(args.output, args.radial, args.angular)
