#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = ["gmsh==4.15.2", "numpy>=2,<3"]
# ///
"""Reproducible public STL inputs for Kutta's external-aero mesher.

Ahmed uses the existing rounded 35-degree body surface generator (no reference
measurement data). FSAE is synthetic, full-scale SI: tapered monocoque and nose,
sidepods, three-element inverted front/rear wings, endplates, floor/diffuser and
four open wheels. Wheels have a 3 mm flattened contact patch at z=0, not a
sub-grid ground clearance. Disjoint named solids are preserved in one ASCII STL.
No volume mesher is used. --kind fsae needs only Python's standard library.
"""
import argparse
import json
import math
from pathlib import Path
import sys


def sub(a, b):
    return tuple(x-y for x, y in zip(a, b))


def cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])


def cap(poly):
    """Ear clipping of a simple 2D cross-section, no convexity assumption."""
    area = sum(a[0]*b[1]-b[0]*a[1] for a, b in zip(poly, poly[1:]+poly[:1]))
    order = list(range(len(poly)))
    if area < 0:
        order.reverse()
    out = []
    orient = lambda a, b, c: (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])
    while len(order) > 3:
        for j, b in enumerate(order):
            a, c = order[j-1], order[(j+1) % len(order)]
            if orient(poly[a], poly[b], poly[c]) <= 1e-14:
                continue
            if any(orient(poly[a], poly[b], poly[p]) >= -1e-14 and
                   orient(poly[b], poly[c], poly[p]) >= -1e-14 and
                   orient(poly[c], poly[a], poly[p]) >= -1e-14
                   for p in order if p not in (a, b, c)):
                continue
            out.append((a, b, c))
            order.pop(j)
            break
        else:
            raise ValueError("non-simple or degenerate section")
    return out + [tuple(order)]


def extrude(poly, ymin, ymax):
    p = [(x, y, z) for y in (ymin, ymax) for x, z in poly]
    n = len(poly)
    faces = []
    for a, b, c in cap(poly):
        faces += [(a, c, b), (a+n, b+n, c+n)]
    for a in range(n):
        b = (a+1) % n
        faces += [(a, b, b+n), (a, b+n, a+n)]
    return p, faces


def box(x0, x1, y0, y1, z0, z1):
    return extrude([(x0,z0),(x1,z0),(x1,z1),(x0,z1)], y0, y1)


def wing(x, z, chord, span, angle):
    # Closed NACA0012 thickness with inverted camber; cosine spacing resolves
    # the nose. Geometric trailing-edge thickness avoids zero-volume exports.
    upper, lower = [], []
    for i in range(33):
        t = (1-math.cos(math.pi*i/32))/2
        h = 5*.12*(.2969*math.sqrt(t)-.126*t-.3516*t*t+.2843*t**3-.1015*t**4)
        camber = -.02*4*t*(1-t)
        for section, sign in ((upper,1),(lower,-1)):
            xx, zz = chord*t, chord*(camber+sign*h)
            section.append((x+xx*math.cos(angle)-zz*math.sin(angle), z+xx*math.sin(angle)+zz*math.cos(angle)))
    return extrude(upper + list(reversed(lower[1:])), -span/2, span/2)


def fsae():
    solids = []
    add = lambda name, geom: solids.append((name, *geom))
    add("monocoque", extrude([(-.65,.14),(.95,.14),(1.18,.4),(.87,.7),(-.25,.63),(-.65,.35)], -.29,.29))
    add("nose", extrude([(-1.25,.2),(-.6,.16),(-.4,.42),(-1.2,.32)], -.13,.13))
    for side in (-1,1):
        a,b = sorted((side*.3,side*.63))
        add("sidepod-"+("left" if side>0 else "right"), box(-.15,.95,a,b,.16,.48))
    add("undertray-diffuser", extrude([(-.75,.065),(.72,.065),(1.36,.23),(1.36,.25),(.72,.085),(-.75,.085)],-.7,.7))
    for which,x,z,chord,span,angle in (("front",-1.5,.12,.32,1.5,7),("rear",.65,.88,.36,1.15,12)):
        for i in range(3):
            add(f"{which}-wing-element-{i+1}",wing(x+i*chord*.81,z+i*(.045 if which=="front" else .075),chord*(1-.12*i),span,math.radians(angle+8*i)))
        for sign in (-1,1):
            y=sign*span/2
            add(f"{which}-wing-endplate-"+("left" if sign>0 else "right"),box(x-.02,x+chord*2.6,y-.004,y+.004,max(.02,z-.08),z+.36))
    for which,x in (("f",-.78),("r",.95)):
        for side,y in (("l",.72),("r",-.72)):
            r=.23
            polygon=[(x+r*math.cos(2*math.pi*i/64),max(0.,r-.003+r*math.sin(2*math.pi*i/64))) for i in range(64)]
            # Remove collinear points on the flattened contact patch.
            polygon=[p for i,p in enumerate(polygon) if not (p[1]==0 and polygon[i-1][1]==0 and polygon[(i+1)%len(polygon)][1]==0)]
            add(f"wheel-{which}{side}",extrude(polygon,y-.105,y+.105))
    return solids


def write(path, solids):
    with path.open("w", newline="\n") as f:
        for name, points, tris in solids:
            f.write(f"solid {name}\n")
            for tri in tris:
                a,b,c = [points[int(i)] for i in tri]
                n = cross(sub(b,a),sub(c,a))
                mag = math.sqrt(sum(x*x for x in n))
                if mag <= 1e-18:
                    raise ValueError(f"degenerate triangle in {name}")
                f.write(" facet normal " + " ".join(f"{x/mag:.12g}" for x in n) + "\n  outer loop\n")
                for p in (a,b,c):
                    f.write("   vertex " + " ".join(f"{v:.12g}" for v in p) + "\n")
                f.write("  endloop\n endfacet\n")
            f.write(f"endsolid {name}\n")


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--kind", choices=["both","ahmed","fsae"],default="both")
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    args.out.mkdir(parents=True,exist_ok=True)
    if args.kind in ("both","ahmed"):
        sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"fsae_aero"))
        import ahmed
        points,tris=ahmed.surface(35,.025)
        write(args.out/"ahmed.stl",[("body",points,tris)])
    if args.kind in ("both","fsae"):
        write(args.out/"fsae.stl",fsae())
    for name,base,minimum,maximum in (("ahmed",.2,.025,.05),("fsae",.5,.04,.1)):
        setup={"v":1,"base_size":base,"surface":{"min_size":minimum,"max_size":maximum},"max_cells":2000000,"layers":{"count":0},"tunnel":{"upstream":3,"downstream":6,"side":2,"top":2,"ground":True}}
        (args.out/f"{name}-coarse.json").write_text(json.dumps(setup,indent=2)+"\n")
    print(json.dumps({"out":str(args.out),"contact_patch":"3 mm flattened wheel bottom at z=0","units":"m"}))


if __name__=="__main__":
    main()
