from __future__ import annotations

import math

from strpg.phmm import PHMMParameters


def test_phmm_transition_and_emission_probabilities_are_normalized() -> None:
    params = PHMMParameters()
    for indel_open, indel_extend in (
        (params.flank_indel_open, params.flank_indel_extend),
        (params.repeat_indel_open, params.repeat_indel_extend),
    ):
        assert math.isclose(
            (1.0 - 2.0 * indel_open) + indel_open + indel_open, 1.0
        )
        assert math.isclose((1.0 - indel_extend) + indel_extend, 1.0)
    assert math.isclose(
        (1.0 - params.mismatch_rate) + 3.0 * params.mismatch_rate / 3.0,
        1.0,
    )
    assert math.isclose(4.0 * params.insertion_background, 1.0)

