from __future__ import annotations

import math

from strpg.models import AlleleRecord
from strpg.phmm import MotifAwarePHMM, PHMMParameters


def test_single_match_uses_match_transition_and_emission() -> None:
    params = PHMMParameters()
    model = MotifAwarePHMM(params)
    candidate = AlleleRecord(
        allele_id="one",
        repeat_count=1,
        sequence="A",
        repeat_sequence="A",
        motif="A",
        repeat_start=0,
        repeat_end=1,
    )
    # Two possible start positions exist, but only start=0 can emit a match.
    match_path = (
        0.5
        * (1.0 - 2.0 * params.repeat_indel_open)
        * (1.0 - params.mismatch_rate)
    )
    insertion_path = (
        0.5 * params.repeat_indel_open * params.insertion_background
    )
    expected = match_path + insertion_path
    assert math.isclose(
        math.exp(model._forward_log_likelihood_python("A", candidate)),
        expected,
        rel_tol=1e-12,
    )
