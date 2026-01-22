"""Entropy utilities for SCARAB."""

from __future__ import annotations

import numpy as np


def renyi_entropy(probabilities, order):
    """Compute the Renyi entropy (base 2) for a probability vector."""
    probs = np.asarray(probabilities, dtype=float).reshape(-1)
    total = probs.sum()
    if total <= 0:
        return float("nan")
    probs = probs / total
    if order < 0:
        raise ValueError("order must be a non-negative real number")
    if order == 0:
        return float(np.log2(np.count_nonzero(probs)))
    if order == 1:
        nonzero = probs > 0
        return float(-np.sum(probs[nonzero] * np.log2(probs[nonzero])))
    if np.isinf(order):
        return float(-np.log2(probs.max()))
    return float(1.0 / (1.0 - order) * np.log2(np.sum(probs**order)))
