from __future__ import annotations

import math

import numpy as np

from strpg.genotyper import _score_genotypes
from strpg.utils import logsumexp
from test_diploid_likelihood_manual_case import alleles


def test_heterozygous_likelihood_is_per_read_probability_mixture() -> None:
    matrix = np.log(np.array([[0.9, 0.1], [0.2, 0.8]]))
    priors = {(0, 0): 1 / 3, (0, 1): 1 / 3, (1, 1): 1 / 3}
    rows = _score_genotypes(matrix, alleles(), priors)
    hetero = next(row for row in rows if (row["i"], row["j"]) == (0, 1))
    expected = sum(
        logsumexp(list(row)) - math.log(2.0) for row in matrix
    )
    assert math.isclose(hetero["log_likelihood"], expected)

