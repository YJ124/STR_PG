from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

try:
    from numba import njit
except ImportError:  # pragma: no cover - exact Python fallback remains available
    njit = None

from .models import AlleleRecord
from .utils import reverse_complement


@dataclass(frozen=True)
class AlignmentSoftmaxParameters:
    """Historical Smith-Waterman likelihood-ablation parameters."""

    match: int = 2
    mismatch: int = -5
    gap: int = -5
    band: int | None = 500
    temperature: float = 12.0

    def validate(self) -> None:
        if self.temperature <= 0:
            raise ValueError("temperature must be positive")
        if self.band is not None and self.band <= 0:
            raise ValueError("band must be positive or None")


def _smith_waterman_score_python(
    query: str,
    target: str,
    match: int = 2,
    mismatch: int = -5,
    gap: int = -5,
    band: int | None = 500,
) -> int:
    """Return the historical local-alignment score without traceback."""
    previous = [0] * (len(target) + 1)
    best = 0
    for i, query_base in enumerate(query, start=1):
        current = [0] * (len(target) + 1)
        low, high = 1, len(target)
        if band is not None:
            low = max(1, i - band)
            high = min(len(target), i + band)
        for j in range(low, high + 1):
            score = max(
                0,
                previous[j - 1] + (match if query_base == target[j - 1] else mismatch),
                previous[j] + gap,
                current[j - 1] + gap,
            )
            current[j] = score
            best = max(best, score)
        previous = current
    return best


def _smith_waterman_score_array(
    query: np.ndarray,
    target: np.ndarray,
    match: int,
    mismatch: int,
    gap: int,
    band: int,
) -> int:
    previous = np.zeros(len(target) + 1, dtype=np.int64)
    best = 0
    for i0 in range(len(query)):
        i = i0 + 1
        current = np.zeros(len(target) + 1, dtype=np.int64)
        low = 1
        high = len(target)
        if band >= 0:
            low = max(1, i - band)
            high = min(len(target), i + band)
        for j in range(low, high + 1):
            substitution = match if query[i0] == target[j - 1] else mismatch
            score = max(
                0,
                previous[j - 1] + substitution,
                previous[j] + gap,
                current[j - 1] + gap,
            )
            current[j] = score
            best = max(best, score)
        previous = current
    return int(best)


_smith_waterman_score_jit = (
    njit(cache=True)(_smith_waterman_score_array) if njit is not None else None
)


def smith_waterman_score(
    query: str,
    target: str,
    match: int = 2,
    mismatch: int = -5,
    gap: int = -5,
    band: int | None = 500,
) -> int:
    """Return the exact historical score, JIT-accelerated when Numba exists."""
    if _smith_waterman_score_jit is None:
        return _smith_waterman_score_python(
            query, target, match, mismatch, gap, band
        )
    query_array = np.frombuffer(query.encode("ascii"), dtype=np.uint8)
    target_array = np.frombuffer(target.encode("ascii"), dtype=np.uint8)
    return int(
        _smith_waterman_score_jit(
            query_array,
            target_array,
            match,
            mismatch,
            gap,
            -1 if band is None else band,
        )
    )


def alignment_softmax_log_likelihood_matrix(
    reads: list[str],
    alleles: list[AlleleRecord],
    params: AlignmentSoftmaxParameters | None = None,
) -> np.ndarray:
    """Return per-read natural-log probabilities from SW-score softmax.

    This ports only the historical score and temperature-softmax likelihood.
    It intentionally omits historical read caps, early stopping, candidate
    pruning, score filters, and no-call rules.
    """
    p = params or AlignmentSoftmaxParameters()
    p.validate()
    matrix = np.empty((len(reads), len(alleles)), dtype=float)
    for read_index, read in enumerate(reads):
        reverse = reverse_complement(read)
        scores = np.asarray(
            [
                max(
                    smith_waterman_score(
                        read, allele.sequence, p.match, p.mismatch, p.gap, p.band
                    ),
                    smith_waterman_score(
                        reverse, allele.sequence, p.match, p.mismatch, p.gap, p.band
                    ),
                )
                for allele in alleles
            ],
            dtype=float,
        )
        shifted = (scores - float(scores.max())) / p.temperature
        log_normalizer = float(np.log(np.exp(shifted).sum()))
        matrix[read_index, :] = shifted - log_normalizer
    return matrix
