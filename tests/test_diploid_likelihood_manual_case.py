from __future__ import annotations

import math

import numpy as np

from strpg.genotyper import _score_genotypes
from strpg.models import AlleleRecord
from strpg.utils import logsumexp


def alleles() -> list[AlleleRecord]:
    return [
        AlleleRecord("a0", 5, "A", "A", "A", 0, 1),
        AlleleRecord("a1", 6, "C", "C", "C", 0, 1),
    ]


def test_diploid_likelihood_manual_case() -> None:
    matrix = np.array([[math.log(0.8), math.log(0.2)]])
    priors = {(0, 0): 1 / 3, (0, 1): 1 / 3, (1, 1): 1 / 3}
    rows = _score_genotypes(matrix, alleles(), priors)
    by_gt = {(row["i"], row["j"]): row for row in rows}
    expected = logsumexp([math.log(0.8), math.log(0.2)]) - math.log(2)
    assert math.isclose(by_gt[(0, 1)]["log_likelihood"], expected)

