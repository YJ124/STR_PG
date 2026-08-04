from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .utils import logsumexp


@dataclass
class MixtureFit:
    weights: np.ndarray
    log_likelihood: float
    iterations: int
    converged: bool


def fit_allele_mixture(
    log_likelihoods: np.ndarray,
    max_iter: int = 200,
    tol: float = 1e-7,
    min_weight: float = 1e-8,
) -> MixtureFit:
    """EM maximum-likelihood fit for per-read allele likelihoods.

    log_likelihoods has shape (n_reads, n_alleles).
    """
    if log_likelihoods.ndim != 2 or log_likelihoods.shape[1] == 0:
        raise ValueError("log_likelihoods must be a non-empty 2D array")
    n_reads, n_alleles = log_likelihoods.shape
    if n_reads == 0:
        return MixtureFit(np.full(n_alleles, 1.0 / n_alleles), -math.inf, 0, True)

    # Deterministic multi-start: uniform plus starts biased to each allele.
    starts = [np.full(n_alleles, 1.0 / n_alleles)]
    for k in range(min(n_alleles, 8)):
        w = np.full(n_alleles, 0.1 / max(1, n_alleles - 1))
        w[k] = 0.9
        w /= w.sum()
        starts.append(w)
    # Include every diploid 0.5/0.5 configuration as an EM start. This
    # guarantees the flexible mixture search begins from the complete
    # diploid model family used by the likelihood-ratio comparison.
    if n_alleles <= 32:
        for i in range(n_alleles):
            for j in range(i, n_alleles):
                w = np.full(n_alleles, min_weight)
                if i == j:
                    w[i] = 1.0
                else:
                    w[i] = 0.5
                    w[j] = 0.5
                w /= w.sum()
                starts.append(w)

    best: MixtureFit | None = None
    for initial in starts:
        weights = initial.copy()
        previous = -math.inf
        converged = False
        for iteration in range(1, max_iter + 1):
            log_joint = log_likelihoods + np.log(np.clip(weights, min_weight, 1.0))[None, :]
            row_norm = np.array([logsumexp(row) for row in log_joint])
            responsibilities = np.exp(log_joint - row_norm[:, None])
            weights = responsibilities.mean(axis=0)
            weights = np.clip(weights, min_weight, None)
            weights /= weights.sum()
            current = float(row_norm.sum())
            if abs(current - previous) <= tol * (1.0 + abs(previous)):
                converged = True
                break
            previous = current
        fit = MixtureFit(weights, current, iteration, converged)
        if best is None or fit.log_likelihood > best.log_likelihood:
            best = fit
    assert best is not None
    return best
