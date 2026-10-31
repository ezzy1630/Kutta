#!/usr/bin/env python3
"""Wrap an Ahmed body mesh (tools/fsae_aero/ahmed.py .msh) in a legacy Fluent .cas.

`kutta solve` needs a Fluent case. This writes the mesh's sections unchanged,
then the V6 mixing elbow's solver settings (SST, Coupled) and one thread
section per Ahmed zone, copied from the elbow zone of the same type. The flow
is not the Ahmed experiment's: it exists so the memory and parity runs go
through the real `kutta solve` path on an external-aero hex mesh.

  python3 tools/memory_scale/ahmed_case.py --mesh runs/mesh/ahmed35_L1.msh \
      --template validation/fluent_examples/mixing_elbow/mixing_elbow.cas.h5 \
      --out runs/mesh/ahmed35_L1.cas
"""
import argparse
import re
import subprocess
import tempfile


def setting(template, name):
    """One /settings string from a CFF case, through h5dump's raw output."""
    with tempfile.NamedTemporaryFile() as raw:
        subprocess.run(["h5dump", "-d", f"/settings/{name}", "-b", "-o", raw.name, template],
                       check=True, stdout=subprocess.DEVNULL)
        return open(raw.name, "rb").read().rstrip(b"\0").decode("ascii")


def sections(text, code):
    """Top-level `(code ...)` sections of an s-expression stream."""
    out, depth, start, i = [], 0, None, 0
    while i < len(text):
        c = text[i]
        if c == '"':
            i = text.index('"', i + 1)
        elif c == "(":
            if depth == 0:
                start = i
            depth += 1
        elif c == ")":
            depth -= 1
            if depth == 0 and re.match(rf"\(\s*{code}\s", text[start:start + 8]):
                out.append(text[start:i + 1])
        i += 1
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mesh", required=True)
    ap.add_argument("--template", required=True, help="CFF case whose settings are copied (V6 elbow)")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    threads = {}
    for t in sections(setting(a.template, "Thread Variables"), 39):
        head = re.match(r"\(39 \((\d+) ([a-z-]+) ([^ ]+) (\d+)\)", t)
        # The first zone of each type is the template for that type.
        threads.setdefault(head.group(2), t[head.end():])
    zones = []
    with open(a.mesh, "rb") as f:
        for line in f:
            m = re.match(rb"\(45 \((\d+) ([a-z-]+) ([^ )]+)\)", line)
            if m:
                zones.append([g.decode() for g in m.groups()])
    version = setting(a.template, "Version")
    with open(a.out, "wb") as out:
        out.write(f'(1 "fluent{version}.0 from tools/memory_scale/ahmed_case.py")\n'.encode())
        # The thread sections below name the zones; the mesh's own (45 ...)
        # names would be duplicates.
        with open(a.mesh, "rb") as f:
            for line in f:
                if not line.startswith(b"(45 "):
                    out.write(line)
        for name in ("Rampant Variables", "Domain Variables"):
            out.write(b"\n" + setting(a.template, name).encode() + b"\n")
        for zid, kind, name in zones:
            out.write(f"(39 ({zid} {kind} {name} 1){threads[kind]}\n".encode())


if __name__ == "__main__":
    main()
