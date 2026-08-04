from __future__ import annotations

import math

import numpy as np

from strpg.genotyper import _score_genotypes
from test_diploid_likelihood_manual_case import alleles


def test_homozygous_likelihood_equals_sum_of_allele_log_likelihoods() -> None:
    matrix = np.log(np.array([[0.8, 0.2], [0.6, 0.4]]))
    priors = {(0, 0): 1 / 3, (0, 1): 1 / 3, (1, 1): 1 / 3}
    rows = _score_genotypes(matrix, alleles(), priors)
    homo = next(row for row in rows if (row["i"], row["j"]) == (0, 0))
    assert math.isclose(homo["log_likelihood"], math.log(0.8) + math.log(0.6))

