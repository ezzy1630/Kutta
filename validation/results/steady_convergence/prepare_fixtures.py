#!/usr/bin/env python3
"""Populate legacy source-root fixture lookups without copying or changing data."""
import hashlib
import json
import os
from pathlib import Path
import subprocess


def main():
    source = Path(os.environ["KUTTA_FLUENT_EXAMPLES"]).resolve(strict=True)
    report = Path(__file__).resolve().parent
    repo = report.parents[2]
    expected = json.loads((report / "metadata.json").read_text())["fluent_fixtures"]["sha256"]
    created, existing = [], []
    for relative, digest in expected.items():
        fixture = source / relative
        if hashlib.sha256(fixture.read_bytes()).hexdigest() != digest:
            raise ValueError(f"External fixture hash mismatch: {fixture}")
        target = repo / "validation/fluent_examples" / relative
        if target.exists():
            if hashlib.sha256(target.read_bytes()).hexdigest() != digest:
                raise ValueError(f"Existing fixture differs; preserving it: {target}")
            existing.append(relative)
            continue
        if target.is_symlink():
            raise ValueError(f"Existing broken fixture link; preserving it: {target}")
        ignored = subprocess.run(["git", "-C", str(repo), "check-ignore", "-q", str(target)])
        if ignored.returncode != 0:
            raise ValueError(f"Refusing to create a non-ignored fixture: {target}")
        target.parent.mkdir(parents=True, exist_ok=True)
        target.symlink_to(fixture)
        created.append(relative)
    print(json.dumps({"source": str(source), "created": created, "existing_identical": existing}, indent=2))


if __name__ == "__main__":
    main()
