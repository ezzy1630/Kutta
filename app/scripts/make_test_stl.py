#!/usr/bin/env python3
"""Write a small, closed, car-like STL for trying New Case from Geometry and
`kutta mesh` without CAD: a body box, front and rear wings and four wheels as
named ASCII solids (one wall zone each), in millimetres, flow along +x, up +z,
wheels touching z = 0.

    make_test_stl.py OUT_DIR [--rear-wing-chord MM] [--split]

--split writes the rear wing to its own file (rear-wing.stl), as a team would
when swapping a part for a design variant; otherwise everything is in car.stl.
"""
import argparse
import math
import os


def box(x0, y0, z0, x1, y1, z1):
    """Twelve outward-facing triangles of an axis-aligned box."""
    p = [(x, y, z) for z in (z0, z1) for y in (y0, y1) for x in (x0, x1)]
    quads = [(0, 2, 3, 1), (4, 5, 7, 6), (0, 1, 5, 4), (2, 6, 7, 3), (0, 4, 6, 2), (1, 3, 7, 5)]
    return [t for a, b, c, d in quads for t in ((p[a], p[b], p[c]), (p[a], p[c], p[d]))]


def wheel(cx, cz, y0, y1, r, n=32):
    """A closed cylinder along y (the axle), centre (cx, cz)."""
    ring = [(cx + r * math.cos(2 * math.pi * k / n), cz + r * math.sin(2 * math.pi * k / n)) for k in range(n)]
    tris = []
    for k in range(n):
        (xa, za), (xb, zb) = ring[k], ring[(k + 1) % n]
        tris += [((xa, y0, za), (xb, y1, zb), (xb, y0, zb)), ((xa, y0, za), (xa, y1, za), (xb, y1, zb))]
        tris += [((cx, y0, cz), (xa, y0, za), (xb, y0, zb)), ((cx, y1, cz), (xb, y1, zb), (xa, y1, za))]
    return tris


def write(path, solids):
    def normal(a, b, c):
        u = [b[i] - a[i] for i in range(3)]
        v = [c[i] - a[i] for i in range(3)]
        n = (u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0])
        m = math.sqrt(sum(x * x for x in n)) or 1.0
        return tuple(x / m for x in n)

    with open(path, "w") as f:
        for name, tris in solids:
            f.write(f"solid {name}\n")
            for a, b, c in tris:
                f.write(" facet normal %g %g %g\n  outer loop\n" % normal(a, b, c))
                for v in (a, b, c):
                    f.write("   vertex %.4f %.4f %.4f\n" % v)
                f.write("  endloop\n endfacet\n")
            f.write(f"endsolid {name}\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out")
    ap.add_argument("--rear-wing-chord", type=float, default=300.0)
    ap.add_argument("--split", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    r = 230.0  # wheel radius
    body = [("body", box(-200, -300, 120, 2400, 300, 650))]
    front = [("front-wing", box(-500, -700, 60, -150, 700, 90))]
    rear = [("rear-wing", box(2300, -450, 950, 2300 + a.rear_wing_chord, 450, 1000))]
    wheels = [("wheel-front", wheel(300, r, 330, 560, r) + wheel(300, r, -560, -330, r)),
              ("wheel-rear", wheel(1830, r, 330, 560, r) + wheel(1830, r, -560, -330, r))]
    if a.split:
        write(os.path.join(a.out, "car.stl"), body + front + wheels)
        write(os.path.join(a.out, "rear-wing.stl"), rear)
    else:
        write(os.path.join(a.out, "car.stl"), body + front + rear + wheels)
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
