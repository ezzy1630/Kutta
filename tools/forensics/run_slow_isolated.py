#!/usr/bin/env python3
"""Run every slow case in a fresh process, preserving its doctest assertions.

This avoids retaining allocator/driver memory across unrelated CFD cases. It
does not change any case, solve budget, reference or gate. The combined XML is
compared by the existing compare_slow_assertions.py, without rewriting INFO or
expression data. macOS physical-footprint monitoring includes shared Metal
memory; output belongs in an ignored build directory.

An interrupted run can resume only from a complete prefix bound to the same
test and CLI binaries. Unreceipted outputs are preserved rather than overwritten.
"""
import argparse
import gc
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import signal
import subprocess
import sys
import time
import xml.etree.ElementTree as ET
from xml.sax.saxutils import quoteattr

from compare_slow_assertions import read_assertions


FOOTPRINT_SOURCE = r"""
#include <libproc.h>
#include <sys/resource.h>
#include <cstdio>
#include <cstdlib>
int main(int argc, char** argv) {
  if (argc != 2) return 1;
  rusage_info_v4 usage{};
  if (proc_pid_rusage(std::atoi(argv[1]), RUSAGE_INFO_V4,
      reinterpret_cast<rusage_info_t*>(&usage))) return 2;
  std::printf("%llu %llu\n", (unsigned long long)usage.ri_phys_footprint,
      (unsigned long long)usage.ri_lifetime_max_phys_footprint);
}
"""


def save(path, value):
    path.write_text(json.dumps(value, indent=2) + "\n")


def stop_group(process):
    if process.poll() is not None:
        return
    try:
        os.killpg(process.pid, signal.SIGTERM)
        process.wait(timeout=5)
    except subprocess.TimeoutExpired:
        os.killpg(process.pid, signal.SIGKILL)
        process.wait()
    except ProcessLookupError:
        process.wait()


def group_memory(helper, group):
    listing = subprocess.check_output(["ps", "-axo", "pid=,pgid="], text=True)
    current, peak = 0, 0
    for line in listing.splitlines():
        pid, pgid = map(int, line.split())
        if pgid != group:
            continue
        result = subprocess.run([str(helper), str(pid)], capture_output=True, text=True)
        if result.returncode == 0:
            live, high = map(int, result.stdout.split())
            current += live
            peak = max(peak, high)
    return current, peak


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, default=Path("build/kutta-tests"))
    parser.add_argument("--cli", type=Path, help="also freeze the CLI used by these tests")
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("--limit-bytes", type=int, default=3_000_000_000)
    parser.add_argument("--resume", action="store_true",
                        help="reuse a verified complete prefix from the same frozen binaries")
    parser.add_argument("--validate-resume", action="store_true",
                        help="validate the saved prefix without starting a solve or changing receipts")
    args = parser.parse_args()
    if sys.platform != "darwin":
        parser.error("physical-footprint monitoring currently requires macOS")
    if args.limit_bytes <= 0:
        parser.error("--limit-bytes must be positive")
    if args.validate_resume and not args.resume:
        parser.error("--validate-resume requires --resume")
    for key in ("KUTTA_FLUENT_EXAMPLES", "KUTTA_RUNS"):
        if not os.environ.get(key):
            parser.error(f"{key} must select the read-only fixture source")

    binary, directory = args.binary.resolve(), args.directory.resolve()
    cli = args.cli.resolve() if args.cli else None
    hashes = {str(p): hashlib.sha256(p.read_bytes()).hexdigest()
              for p in (binary, cli) if p is not None}
    shards = directory / "shards"
    ledger = directory / "isolated-results.json"
    if args.resume:
        if not directory.is_dir() or not shards.is_dir() or not ledger.is_file():
            raise ValueError("resume requires an existing shard directory and receipt")
    else:
        directory.mkdir(parents=True, exist_ok=True)
        shards.mkdir()
    receipt = json.loads(ledger.read_text()) if args.resume else None
    if receipt is not None:
        if (receipt.get("complete") is not False or receipt.get("frozen_binaries") != hashes
                or receipt.get("binary_sha256") != hashes[str(binary)]
                or receipt.get("limit_bytes") != args.limit_bytes):
            raise ValueError("resume requires an incomplete receipt from these exact binaries and resource limit")
    env = dict(os.environ, MallocSpaceEfficient="1")
    env.pop("KUTTA_MEMORY_PHASES", None)  # CLI INFO is part of assertion identity.
    manifest = directory / "case-manifest.xml"
    listing = subprocess.run([str(binary), "-ts=slow", "-ltc", "-r=xml"],
                             stdout=subprocess.PIPE, env=env, check=True)
    if not args.resume:
        manifest.write_bytes(listing.stdout)
    manifest_root = ET.parse(manifest).getroot()
    version = manifest_root.get("version", "unknown")
    names = [c.get("name") for c in manifest_root.iter("TestCase")
             if c.get("skipped") != "true"]
    if not names or len(names) != len(set(names)):
        raise ValueError("slow-case manifest is empty or contains duplicate names")
    listed_names = [c.get("name") for c in ET.fromstring(listing.stdout).iter("TestCase")
                    if c.get("skipped") != "true"]
    if listed_names != names:
        raise ValueError("current binary's slow-case manifest differs from the saved manifest")

    totals = {"successes": 0, "failures": 0}
    completed_prefix = 0
    if receipt is not None:
        if receipt.get("case_count") != len(names):
            raise ValueError("saved receipt has a different slow-case count")
        for row in receipt["cases"]:
            index = completed_prefix + 1
            if (index > len(names) or row.get("index") != index or row.get("case") != names[index - 1]
                    or row.get("complete") is not True or row.get("exit") not in (0, 1)
                    or row.get("resource_limit_exceeded") is not False
                    or row.get("reported_peak_bytes") is None
                    or max(row["sampled_peak_bytes"], row["reported_peak_bytes"]) > args.limit_bytes):
                raise ValueError(f"case {index} is not a verified complete prefix entry")
            expected_command = ["/usr/bin/time", "-l", str(binary), "-ts=slow", "-r=xml",
                                "--success", f"--first={index}", f"--last={index}"]
            if row.get("command") != expected_command:
                raise ValueError(f"case {index} was run with different command options")
            part = shards / f"case-{index:03d}.xml"
            fragment = part.with_suffix(".fragment.xml")
            actual = read_assertions(fragment)
            if list(actual) != [names[index - 1]]:
                raise ValueError(f"saved fragment {index} has a different case identity")
            del actual
            selected = [c for c in ET.parse(part).getroot().iter("TestCase")
                        if c.get("skipped") != "true"]
            if len(selected) != 1 or ET.tostring(selected[0], encoding="utf-8") != fragment.read_bytes():
                raise ValueError(f"saved fragment {index} differs from its original test output")
            summary = selected[0].find("OverallResultsAsserts")
            if dict(summary.attrib) != row.get("assertions"):
                raise ValueError(f"saved assertion totals differ for case {index}")
            for key in totals:
                totals[key] += int(summary.get(key, 0))
            completed_prefix = index
            del selected, summary
            gc.collect()
    if args.validate_resume:
        print(f"VERIFIED complete prefix {completed_prefix}/{len(names)}; no solve started", flush=True)
        return 0

    helper_source, helper = directory / "footprint-pid.cpp", directory / "footprint-pid"
    if not args.resume or not helper.exists():
        helper_source.write_text(FOOTPRINT_SOURCE)
        try:
            subprocess.run([*shlex.split(os.environ.get("CXX", "c++")), "-O2",
                            str(helper_source), "-o", str(helper)], check=True)
        finally:
            helper_source.unlink()
    if receipt is None:
        receipt = {"mode": "one fresh process per slow test case",
                   "binary_sha256": hashes[str(binary)], "frozen_binaries": hashes,
                   "case_count": len(names), "limit_bytes": args.limit_bytes,
                   "complete": False, "cases": []}
    else:
        receipt.setdefault("resumptions", []).append({
            "utc_epoch": time.time(), "verified_complete_prefix": completed_prefix,
            "prior_receipt_sha256": hashlib.sha256(ledger.read_bytes()).hexdigest(),
            "runner_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        })
    save(ledger, receipt)

    def interrupted(signum, frame):
        raise KeyboardInterrupt(f"interrupted by signal {signum}")
    signal.signal(signal.SIGTERM, interrupted)
    for index, name in enumerate(names, 1):
        if index <= completed_prefix:
            continue
        for path, expected in hashes.items():
            if hashlib.sha256(Path(path).read_bytes()).hexdigest() != expected:
                raise ValueError(f"frozen binary changed during the run: {path}")
        part = shards / f"case-{index:03d}.xml"
        errors = part.with_suffix(".stderr")
        if any(p.exists() for p in (part, errors, part.with_suffix(".fragment.xml"))):
            raise ValueError(f"unreceipted shard {index} exists; preserve it before explicitly retrying")
        command = ["/usr/bin/time", "-l", str(binary), "-ts=slow", "-r=xml",
                   "--success", f"--first={index}", f"--last={index}"]
        print(f"START {index}/{len(names)} {name}", flush=True)
        started = time.monotonic()
        maximum, exceeded = 0, False
        with part.open("wb") as output, errors.open("wb") as stderr:
            process = subprocess.Popen(command, env=env, stdout=output, stderr=stderr,
                                       start_new_session=True)
            try:
                heartbeat = started
                while process.poll() is None:
                    live, high = group_memory(helper, process.pid)
                    maximum = max(maximum, live, high)
                    if maximum > args.limit_bytes:
                        exceeded = True
                        stop_group(process)
                        break
                    if time.monotonic() - heartbeat >= 60:
                        print(f"RUNNING {index}/{len(names)} peak_bytes={maximum}", flush=True)
                        heartbeat = time.monotonic()
                    time.sleep(1)
                code = process.wait()
            finally:
                stop_group(process)
        peaks = re.findall(r"(?m)^\s*(\d+)\s+peak memory footprint", errors.read_text())
        reported_peak = int(peaks[-1]) if peaks else None
        if reported_peak is not None and reported_peak > args.limit_bytes:
            exceeded = True
        row = {"index": index, "case": name, "command": command,
               "exit": code, "sampled_peak_bytes": maximum,
               "reported_peak_bytes": reported_peak,
               "resource_limit_exceeded": exceeded,
               "elapsed_seconds": time.monotonic() - started}
        receipt["cases"].append(row)
        save(ledger, receipt)
        if exceeded:
            print(f"STOP resource limit: {name}, {max(maximum, reported_peak or 0)} bytes", flush=True)
            return 3
        actual = read_assertions(part)  # Use the unchanged canonical parser.
        if list(actual) != [name]:
            raise ValueError(f"shard {index} executed {list(actual)}, expected {name}")
        del actual
        tree = ET.parse(part)
        cases = [c for c in tree.getroot().iter("TestCase") if c.get("skipped") != "true"]
        if len(cases) != 1 or cases[0].get("name") != name:
            raise ValueError(f"shard {index} case identity changed")
        summary = cases[0].find("OverallResultsAsserts")
        for key in totals:
            totals[key] += int(summary.get(key, 0))
        fragment = part.with_suffix(".fragment.xml")
        ET.ElementTree(cases[0]).write(fragment, encoding="utf-8")
        row["assertions"] = dict(summary.attrib)
        row["complete"] = True
        save(ledger, receipt)
        del summary, cases, tree
        gc.collect()
        print(f"DONE {index}/{len(names)} exit={code} peak_bytes={maximum}", flush=True)

    combined = directory / "slow.xml"
    with combined.open("wb") as output:
        output.write(b'<?xml version="1.0" encoding="UTF-8"?>\n')
        output.write(f'<doctest binary={quoteattr(str(binary))} version={quoteattr(version)}>\n'.encode())
        output.write(f'<Execution mode="isolated-per-test-case" cases="{len(names)}"/>\n'.encode())
        output.write(b'<TestSuite name="slow">\n')
        for index in range(1, len(names) + 1):
            with (shards / f"case-{index:03d}.fragment.xml").open("rb") as fragment:
                shutil.copyfileobj(fragment, output)
        output.write(b'\n</TestSuite>\n')
        output.write(ET.tostring(ET.Element("OverallResultsAsserts", {k: str(v) for k, v in totals.items()})))
        output.write(b'\n</doctest>\n')
    actual = read_assertions(combined)
    if list(actual) != names:
        raise ValueError("combined XML did not preserve the complete case manifest")
    receipt.update(complete=True, assertions=totals,
                   combined_sha256=hashlib.sha256(combined.read_bytes()).hexdigest())
    save(ledger, receipt)
    print(f"COMPLETE {len(names)} cases, {totals}", flush=True)
    return 1 if totals["failures"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
