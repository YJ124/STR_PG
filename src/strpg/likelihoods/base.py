"""Backend contract for read-by-allele likelihood calculation."""

from __future__ import annotations

from abc import ABC, abstractmethod

import numpy as np

from ..models import AlleleRecord, LocusRecord


class AlleleLikelihoodBackend(ABC):
    """Score a fixed ordered read set against a fixed ordered allele set.

    Implementations return natural-log-compatible values where larger values
    indicate stronger support. They must not select reads, alter candidates,
    add priors, enumerate genotypes, or inspect truth.
    """

    name: str

    @abstractmethod
    def score_reads(
        self,
        reads: list[str],
        candidate_alleles: list[AlleleRecord],
        locus: LocusRecord | None = None,
    ) -> np.ndarray:
        """Return a `(len(reads), len(candidate_alleles))` matrix."""
