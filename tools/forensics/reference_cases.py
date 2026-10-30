#!/usr/bin/env python3
"""Restore, verify and inspect the licensed Fluent research corpus. Parse only.

No downloaded journal, UDF, Python, Scheme or project command is executed.
A successful one-step probe tests adapter admission, never convergence.
"""
import argparse
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import subprocess
import urllib.request
import zipfile

ROOT = Path(__file__).resolve().parents[2]
CORPUS = ROOT / "validation/fluent_examples"
MANIFEST = ROOT / "validation/fluent_reference_manifest.tsv"


def ledger():
    with MANIFEST.open() as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    seen = set()
    for row in rows:
        path = Path(row["path"])
        if path.is_absolute() or ".." in path.parts or row["path"] in seen:
            raise ValueError("invalid or duplicate ledger path")
        seen.add(row["path"])
        row["bytes"] = int(row["bytes"])
    return rows


def valid(row, data):
    return len(data) == row["bytes"] and hashlib.sha256(data).hexdigest() == row["sha256"]


def restore(row):
    # Only pre-reviewed manifest URLs are used. Extract one named member rather
    # than unpacking arbitrary paths, permissions, symlinks or executables.
    with urllib.request.urlopen(row["url"], timeout=120) as response:
        raw = response.read(256 * 1024 * 1024 + 1)
    if len(raw) > 256 * 1024 * 1024:
        raise ValueError("download exceeds corpus size bound")
    if row["transform"] == "direct":
        data = raw
    elif row["transform"] == "gzip":
        with gzip.GzipFile(fileobj=io.BytesIO(raw)) as stream:
            data = stream.read(row["bytes"] + 1)
    elif row["transform"] == "zip":
        with zipfile.ZipFile(io.BytesIO(raw)) as archive:
            info = archive.getinfo(row["member"])
            if info.file_size != row["bytes"]:
                raise ValueError("archive member size mismatch")
            data = archive.read(info)
    else:
        raise ValueError("unknown corpus transform")
    if not valid(row, data):
        raise ValueError("download hash or byte count mismatch: " + row["path"])
    path = CORPUS / row["path"]
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".download-tmp")
    temporary.write_bytes(data)
    temporary.replace(path)


def pairs(case):
    data = case.with_name(case.name.replace(".cas", ".dat"))
    if data.exists():
        return [data]
    if case.name == "elbow-2.cas.h5":
        return sorted(case.parent.glob("elbow-2-*.dat.h5"))
    return []


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--download", action="store_true", help="restore missing files; reject changed existing files")
    parser.add_argument("--inspect", action="store_true")
    parser.add_argument("--probe", action="store_true", help="strict then named WFGC one-step solve admission")
    parser.add_argument("--kutta", type=Path, default=ROOT / "build/kutta")
    parser.add_argument("--output", type=Path, default=ROOT / "build/fluent-reference-audit")
    args = parser.parse_args()
    rows = ledger()
    missing = []
    for row in rows:
        path = CORPUS / row["path"]
        if not path.exists() and args.download:
            restore(row)
        if not path.exists():
            missing.append(row["path"])
        elif not valid(row, path.read_bytes()):
            raise ValueError("existing file differs from ledger: " + row["path"])
    print(f"Verified {len(rows)-len(missing)}/{len(rows)} files; absent {len(missing)}")
    if not args.inspect and not args.probe:
        return
    args.output.mkdir(parents=True, exist_ok=True)
    results = []
    for row in rows:
        case = CORPUS / row["path"]
        if not case.exists() or not case.name.endswith((".cas", ".cas.h5")):
            continue
        saved = pairs(case)
        for data in saved or [None]:
            record = {"case": row["path"], "data": str(data.relative_to(CORPUS)) if data else None}
            if args.inspect:
                cmd = [str(args.kutta), "inspect", str(case)]
                if data:
                    cmd += ["--data", str(data)]
                result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
                record["inspect"] = {"exit": result.returncode, "output": result.stdout + result.stderr}
            if args.probe and data:
                record["probes"] = []
                for approximation in [None, "wfgc"]:
                    destination = args.output / "probes" / row["path"] / data.name / str(approximation)
                    cmd = [str(args.kutta), "solve", str(case), "--data", str(data), "--backend", "cpu", "--iterations", "1", "--out", str(destination)]
                    if approximation:
                        cmd += ["--approximate", approximation]
                    result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
                    record["probes"].append({"approximation": approximation, "exit": result.returncode, "output": result.stdout + result.stderr})
                    if result.returncode == 0:
                        break
            results.append(record)
            print(row["path"], record["data"], flush=True)
    (args.output / "audit.json").write_text(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
