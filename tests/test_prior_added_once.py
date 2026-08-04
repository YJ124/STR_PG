from __future__ import annotations

import math

import numpy as np

from strpg.genotyper import _score_genotypes
from test_diploid_likelihood_manual_case import alleles


def test_prior_is_added_once_after_all_reads() -> None:
    matrix = np.log(np.array([[0.8, 0.2], [0.7, 0.3], [0.9, 0.1]]))
    priors = {(0, 0): 0.5, (0, 1): 0.3, (1, 1): 0.2}
    rows = _score_genotypes(matrix, alleles(), priors)
    homo = next(row for row in rows if (row["i"], row["j"]) == (0, 0))
    expected_likelihood = sum(math.log(value) for value in (0.8, 0.7, 0.9))
    assert math.isclose(homo["log_likelihood"], expected_likelihood)
    assert math.isclose(
        homo["log_score"], expected_likelihood + math.log(0.5)
    )

