from __future__ import annotations

import math

from strpg.models import AlleleRecord
from strpg.phmm import MotifAwarePHMM
from test_phmm_forward_matches_bruteforce import brute_force_probability


def test_identical_flank_base_has_equal_free_start_weight() -> None:
    model = MotifAwarePHMM()
    candidate = AlleleRecord(
        allele_id="boundary",
        repeat_count=0,
        sequence="ACA",
        repeat_sequence="",
        motif="C",
        repeat_start=1,
        repeat_end=1,
    )
    # The brute-force enumerator includes both exchangeable A starts plus all
    # insertion/deletion paths and does not average mutually exclusive ends.
    observed = math.exp(
        model._forward_log_likelihood_python("A", candidate)
    )
    expected = brute_force_probability(
        "A", candidate, model.params
    )
    assert math.isclose(observed, expected, rel_tol=1e-12)
