from __future__ import annotations

import numpy as np

from strpg.genotyper import _score_genotypes
from test_diploid_likelihood_manual_case import alleles


def test_candidate_matrix_columns_remain_aligned_with_allele_indices() -> None:
    candidates = alleles()
    matrix = np.log(np.array([[0.01, 0.99], [0.02, 0.98]]))
    priors = {(0, 0): 1 / 3, (0, 1): 1 / 3, (1, 1): 1 / 3}
    best = _score_genotypes(matrix, candidates, priors)[0]
    assert candidates[best["i"]].allele_id == "a1"
    assert candidates[best["j"]].allele_id == "a1"
