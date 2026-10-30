#!/usr/bin/env python3
"""Expand bundled TMR grids for the ASCII-only Plot3D reader; no downloads."""
import gzip
import hashlib
from pathlib import Path

root = Path(__file__).resolve().parents[1] / "validation/data/nasa_tmr"
# SOURCE.md pins compressed inputs. Verify them before creating derived files.
grids = {
    "2DBC_bump/grids/bump_2levelsdown_353x161.p2dfmt.gz":
        "e0b49f79da69ff2cdab2f0d7f8628af2a7823b6afb0e1d48d587af0186ba6689",
    "2DZP_flatplate/grids/flatplate_clust2_1leveldown_273x193.p2dfmt.gz":
        "40be4c0a92b75892189263a5c3efb63fa112030912480427e1e0787c6cb97955",
}
for relative, expected in grids.items():
    source = root / relative
    data = source.read_bytes()
    if hashlib.sha256(data).hexdigest() != expected:
        raise SystemExit(f"TMR checksum mismatch: {source}")
    target = source.with_suffix("")
    expanded = gzip.decompress(data)
    if target.exists() and target.read_bytes() != expanded:
        raise SystemExit(f"Refusing to overwrite different grid: {target}")
    target.write_bytes(expanded)
    print(f"{target}: {len(expanded)} bytes")
