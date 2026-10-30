#!/usr/bin/env python3
"""Parse stored mesh types, WFGC flags and residuals; requires h5py/numpy.

Stored Scheme stays text. No downloaded code is evaluated. Numerical results
come from kutta inspect/solve, not this metadata reader.
"""
from collections import Counter
import json
from pathlib import Path
import re
import struct

import h5py

ROOT = Path(__file__).resolve().parents[2]
FLAGS = ["recon/cell-lsf?", "recon/wfgc-jump-enhanced?", "recon/node-lsq?",
         "rp-unsteady?", "rp-lam?", "number-of-iterations"]


def flags(text):
    return {key: re.findall(r"\(" + re.escape(key) + r"\s+([^\s()]+)\)", text)
            for key in FLAGS}


def hdf5_metadata(path):
    result = {}
    with h5py.File(path) as source:
        settings = {}
        if "settings" in source:
            for key, value in source["settings"].items():
                if isinstance(value, h5py.Dataset):
                    settings[key] = value[()].tobytes().decode(errors="replace").rstrip("\0")
        result["flags"] = flags(settings.get("Rampant Variables", ""))
        result["settings"] = {key: settings[key] for key in
                              ["Origin", "Solver", "Version", "Case File"] if key in settings}
        counts = Counter()
        if "meshes/1/cells/ctype" in source:
            for section in source["meshes/1/cells/ctype"].values():
                if "cell-types" in section:
                    counts.update(map(int, section["cell-types"][()].ravel()))
                elif "elementType" in section.attrs:
                    kind = int(section.attrs["elementType"][0])
                    counts[kind] += int(section.attrs["maxId"][0] - section.attrs["minId"][0] + 1)
        result["cell_types"] = dict(counts)
        if "meshes/1" in source:
            result["mesh_attributes"] = {
                key: source["meshes/1"].attrs[key].tolist() for key in
                ["cellCount", "faceCount", "nodeCount", "dimension"]
                if key in source["meshes/1"].attrs}
        if "results/residuals" in source:
            result["residuals"] = {}
            for equation, group in source["results/residuals"].items():
                if "data" in group and len(group["data"]):
                    result["residuals"][equation] = {
                        "last_iteration": float(group["iterations"][-1]),
                        "last": group["data"][-1].tolist()}
        for location in ["faces", "cells"]:
            group = "results/1/phase-1/" + location
            if group in source:
                result[location + "_fields"] = list(source[group])
    return result


def legacy_metadata(path):
    data = path.read_bytes()
    text = data.decode("latin1")
    declarations = re.findall(r"\((?:12|2012|3012)\s*\(([0-9a-fA-F ]+)\)", text)
    counts = Counter()
    for declaration in declarations:
        header = [int(token, 16) for token in declaration.split()]
        if len(header) == 5 and header[0] != 0 and header[4] != 0:
            counts[header[4]] += header[2] - header[1] + 1
    # Mixed legacy types are 32-bit ints even in double-precision section 3012.
    # Parse only numeric headers; the file length bounds the binary read.
    for match in re.finditer(rb"\((2012|3012)\s*\(([0-9a-fA-F ]+)\)\s*\(", data):
        header = [int(token, 16) for token in match[2].split()]
        if len(header) != 5 or header[0] == 0 or header[4] != 0:
            continue
        count = header[2] - header[1] + 1
        start = match.end()
        if data[start:start+2] == b"\r\n":
            start += 2
        elif data[start:start+1] in [b"\n", b" "]:
            start += 1
        if count < 0 or 4 * count > len(data) - start:
            raise ValueError("legacy mixed type range exceeds file")
        counts.update(struct.unpack_from("<" + str(count) + "i", data, start))
    return {"flags": flags(text), "cell_types": dict(counts),
            "legacy_cell_declarations": declarations,
            "legacy_data_headers": re.findall(r"\((?:300|2300|3300)\s*\(([^)]*)\)", text)[:12]}


def main():
    rows = []
    for path in sorted((ROOT / "validation/fluent_examples").rglob("*")):
        if not path.name.endswith((".cas", ".dat", ".cas.h5", ".dat.h5")):
            continue
        row = {"file": str(path.relative_to(ROOT))}
        try:
            row.update(hdf5_metadata(path) if h5py.is_hdf5(path) else legacy_metadata(path))
        except (OSError, ValueError, KeyError, struct.error) as error:
            row["error"] = str(error)
        rows.append(row)
    output = ROOT / "validation/results/fluent_reference_cases/metadata.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(rows, indent=2))
    print(f"Parsed {len(rows)} files: {output}")


if __name__ == "__main__":
    main()
