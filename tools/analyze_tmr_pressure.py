#!/usr/bin/env python3
"""Independently inspect captured/replayed V5 pressure matrices (requires numpy)."""
import argparse
import gzip
import math
from pathlib import Path
import numpy as np


def load(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as stream:
        cells, iteration = map(int, next(stream).split())
        matrix = np.zeros((cells, cells))
        for line in stream:
            row, column, value = line.split()
            matrix[int(row), int(column)] = float(value)
    return matrix, iteration


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("path", type=Path, nargs="?", default=Path(__file__).resolve().parents[1] / "validation/results/tmr_v5/bump_pressure_iteration18.coo.gz")
parser.add_argument("--compare", type=Path, help="verify replay retains the original off-diagonals")
parser.add_argument("--require-spd", action="store_true", help="fail unless symmetric, diagonally dominant and Cholesky-positive")
args = parser.parse_args()
matrix, iteration = load(args.path)
cells = len(matrix)
diagonal = np.diag(matrix)
surplus = np.array([math.fsum(row[row != 0]) for row in matrix])
symmetry = np.max(np.abs(matrix - matrix.T))
print(f"completed iterations={iteration}, cells={cells}")
print(f"symmetry error={symmetry:.12g}")
print(f"diagonal min/max={diagonal.min():.12g}/{diagonal.max():.12g}")
print(f"negative row sums={np.count_nonzero(surplus < 0)}, minimum row surplus={surplus.min():.12g}")
eigenvalues = np.linalg.eigvalsh(matrix)
print(f"smallest eigenvalues={eigenvalues[:5]}")
if eigenvalues[0] > 0:
    print(f"unscaled condition estimate={eigenvalues[-1] / eigenvalues[0]:.12g}")
positive = False
try:
    np.linalg.cholesky(matrix)
    positive = True
    print("Cholesky: PASS")
except np.linalg.LinAlgError:
    print("Cholesky: FAIL (stored matrix is not positive definite)")
if np.all(diagonal > 0):
    scale = np.sqrt(diagonal)
    scaled = matrix / scale[:, None] / scale[None, :]
    eigenvalues = np.linalg.eigvalsh(scaled)
    print(f"diagonal-congruence eigenvalue min/max={eigenvalues[0]:.12g}/{eigenvalues[-1]:.12g}")
    if eigenvalues[0] > 0:
        print(f"diagonal-congruence condition estimate={eigenvalues[-1] / eigenvalues[0]:.12g}")
    else:
        print("diagonal congruence retains an indefinite operator")
same_edges = True
if args.compare:
    original, _ = load(args.compare)
    difference = matrix - original
    change = np.diag(difference).copy()
    np.fill_diagonal(difference, 0)
    same_edges = not np.any(difference)
    print(f"off-diagonals unchanged={same_edges}, changed diagonals={np.count_nonzero(change)}")
    print(f"diagonal change min/max={change.min():.12g}/{change.max():.12g}")
if args.require_spd and not (positive and symmetry == 0 and np.all(surplus >= 0) and same_edges):
    raise SystemExit(1)
