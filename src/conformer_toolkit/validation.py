"""Shared validation for conformer data."""

from __future__ import annotations

import numpy as np


def validate_rmsd_matrix(rmsd_matrix: np.ndarray) -> np.ndarray:
    """Return an RMSD matrix as float64 after validating its shape and values."""

    matrix = np.asarray(rmsd_matrix, dtype=np.float64)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("rmsd_matrix must be a square N x N array")
    if matrix.shape[0] < 2:
        raise ValueError("rmsd_matrix must contain at least two conformers")
    if not np.all(np.isfinite(matrix)):
        raise ValueError("rmsd_matrix must contain only finite values")
    if np.any(matrix < 0.0):
        raise ValueError("rmsd_matrix cannot contain negative distances")
    return matrix

