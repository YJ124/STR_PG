"""Optional experimental motif-aware profile-HMM likelihood backend."""

from __future__ import annotations

import numpy as np

from ..models import AlleleRecord, LocusRecord
from ..phmm import MotifAwarePHMM, PHMMParameters
from .base import AlleleLikelihoodBackend


class PHMMLikelihoodBackend(AlleleLikelihoodBackend):
    """Expose the audited pHMM through the common backend contract."""

    name = "phmm"

    def __init__(self, params: PHMMParameters | None = None):
        self.model = MotifAwarePHMM(params)

    def score_reads(
        self,
        reads: list[str],
        candidate_alleles: list[AlleleRecord],
        locus: LocusRecord | None = None,
    ) -> np.ndarray:
        matrix = np.empty(
            (len(reads), len(candidate_alleles)), dtype=float
        )
        for read_index, read in enumerate(reads):
            for allele_index, allele in enumerate(candidate_alleles):
                matrix[read_index, allele_index] = (
                    self.model.orientation_marginal_log_likelihood(read, allele)
                )
        return matrix
