"""Make the explicitly WFGC-disabled local CLI diagnostic; requires h5py."""
import hashlib
import os
from pathlib import Path
import shutil

import h5py

expected_source = "8145813f50f989a521ed8d1b02c891c7ca22ccc5fafa35a4eeba011949f19cbc"
expected_target = "b3b550e8b9c2a2902c45319d7543cc5eada682367b3e53fb8427a2748e69a0ff"
source = Path(os.environ["KUTTA_FLUENT_EXAMPLES"]) / "mixing_elbow/mixing_elbow.cas.h5"
target = Path("build/elbow-hybrid-cold.cas.h5")
digest = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
assert digest(source) == expected_source
assert source.resolve() != target.resolve()
assert not target.exists() or digest(target) in {expected_source, expected_target}
target.parent.mkdir(parents=True, exist_ok=True)
shutil.copyfile(source, target)
with h5py.File(target, "r+") as case:
    variables = case["settings/Rampant Variables"]
    raw = variables[0]
    assert raw.count(b"(recon/cell-lsf? #t)") == 1
    assert b"(hyb-init? #t)" in raw
    variables[0] = raw.replace(b"(recon/cell-lsf? #t)", b"(recon/cell-lsf? #f)")
assert digest(source) == expected_source
assert digest(target) == expected_target
print(f"Local WFGC-disabled diagnostic, stored hybrid controls retained: {target}")
print(f"SHA256 {expected_target}")
