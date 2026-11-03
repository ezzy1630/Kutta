#!/usr/bin/env python3
"""Turn bench_memory KUTTA_MEMORY_BUFFERS output into one CSV row per buffer.

Use --revision for the source that built the binary; line numbers are allocation
provenance, not a guess based on matching a buffer's length to a mesh count.
The snapshot includes each live solver-owned allocation exactly once. Backend
reduction/command storage outside Device::alloc is excluded (see the report).
"""
import argparse
import collections
import csv
import pathlib
import shlex
import subprocess
import sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("log", type=pathlib.Path)
    ap.add_argument("--revision", help="git revision of the allocating source; default is working tree")
    ap.add_argument("--phase", default="iterations")
    a = ap.parse_args()
    lines = a.log.read_text().splitlines()
    total = next((dict(word.split("=", 1) for word in line.split()[1:])
                  for line in reversed(lines) if line.startswith("MEMORY_TOTAL ")), None)
    if total is None:
        raise SystemExit("missing MEMORY_TOTAL: run did not finish")
    cells, faces = int(total["cells"]), int(total["faces"])
    records = []
    sources = {}
    for line in lines:
        if not line.startswith("MEMORY_BUFFER "):
            continue
        row = dict(word.split("=", 1) for word in shlex.split(line)[1:])
        if row["phase"] != a.phase:
            continue
        file = row["file"]
        if file not in sources:
            text = (subprocess.check_output(["git", "show", f"{a.revision}:{file}"], text=True)
                    if a.revision else pathlib.Path(file).read_text())
            sources[file] = text.splitlines()
        row["allocation"] = sources[file][int(row["line"]) - 1].strip()
        records.append(row)
    if not records:
        raise SystemExit("no matching buffer snapshot")
    records.sort(key=lambda r: (r["file"], int(r["line"]), int(r["column"]), int(r["bytes"])))
    writer = csv.writer(sys.stdout)
    writer.writerow(["buffer", "source", "column", "bytes", "bytes_per_cell", "bytes_per_face", "allocation"])
    ids = collections.Counter()
    for row in records:
        source = f'{row["file"]}:{row["line"]}'
        ids[source] += 1
        writer.writerow([f'{source}#{ids[source]}', source, row["column"], row["bytes"],
                         format(int(row["bytes"])/cells,".12g"), format(int(row["bytes"])/faces,".12g"), row["allocation"]])


if __name__ == "__main__":
    main()
