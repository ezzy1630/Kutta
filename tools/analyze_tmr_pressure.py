#!/usr/bin/env python3
"""Independently inspect the captured V5 pressure matrix (requires numpy)."""
import gzip
from pathlib import Path
import numpy as np

path = Path(__file__).resolve().parents[1] / "validation/results/tmr_v5/bump_pressure_iteration18.coo.gz"
with gzip.open(path, "rt") as stream:
    cells, iteration = map(int, next(stream).split())
    matrix = np.zeros((cells, cells))
    for line in stream:
        row, column, value = line.split()
        matrix[int(row), int(column)] = float(value)
diagonal = np.diag(matrix)
print(f"completed iterations={iteration}, cells={cells}")
print(f"symmetry error={np.max(np.abs(matrix - matrix.T)):.12g}")
print(f"diagonal min/max={diagonal.min():.12g}/{diagonal.max():.12g}")
print(f"negative row sums={np.count_nonzero(matrix.sum(axis=1) < 0)}")
print(f"smallest eigenvalues={np.linalg.eigvalsh(matrix)[:5]}")
try:
    np.linalg.cholesky(matrix)
    print("Cholesky: PASS")
except np.linalg.LinAlgError:
    print("Cholesky: FAIL (stored matrix is not positive definite)")
