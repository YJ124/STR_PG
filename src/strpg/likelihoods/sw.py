"""Production Smith-Waterman allele-likelihood backend."""

from __future__ import annotations

import numpy as np

from ..alignment import (
    AlignmentSoftmaxParameters,
    alignment_softmax_log_likelihood_matrix,
)
from ..models import AlleleRecord, LocusRecord
from .base import AlleleLikelihoodBackend


class SWLikelihoodBackend(AlleleLikelihoodBackend):
    """Convert bidirectional local-alignment scores to log probabilities."""

    name = "sw"

    def __init__(self, params: AlignmentSoftmaxParameters | None = None):
        self.params = params or AlignmentSoftmaxParameters()
        self.params.validate()

    def score_reads(
        self,
        reads: list[str],
        candidate_alleles: list[AlleleRecord],
        locus: LocusRecord | None = None,
    ) -> np.ndarray:
        return alignment_softmax_log_likelihood_matrix(
            reads, candidate_alleles, self.params
        )
