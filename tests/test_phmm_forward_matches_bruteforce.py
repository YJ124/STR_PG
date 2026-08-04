from __future__ import annotations

import math
from functools import lru_cache

from strpg.models import AlleleRecord
from strpg.phmm import MotifAwarePHMM, PHMMParameters


def allele(sequence: str, repeat_start: int = 0, repeat_end: int = 0) -> AlleleRecord:
    return AlleleRecord(
        allele_id="toy",
        repeat_count=0,
        sequence=sequence,
        repeat_sequence=sequence[repeat_start:repeat_end],
        motif="A",
        repeat_start=repeat_start,
        repeat_end=repeat_end,
    )


def brute_force_probability(
    read: str, candidate: AlleleRecord, params: PHMMParameters
) -> float:
    """Enumerate every finite M/I/D state path for a tiny semi-global case.

    Start positions are mutually exclusive and have a uniform prior. End
    positions are mutually exclusive outcomes and are summed, not averaged.
    """

    template = candidate.sequence
    n, m = len(read), len(template)

    def transition(position: int) -> tuple[float, float, float, float, float]:
        in_repeat = candidate.repeat_start <= position < candidate.repeat_end
        open_prob = (
            params.repeat_indel_open if in_repeat else params.flank_indel_open
        )
        extend_prob = (
            params.repeat_indel_extend if in_repeat else params.flank_indel_extend
        )
        return (
            1.0 - 2.0 * open_prob,
            open_prob,
            open_prob,
            1.0 - extend_prob,
            extend_prob,
        )

    def match_emission(read_base: str, template_base: str) -> float:
        return (
            1.0 - params.mismatch_rate
            if read_base == template_base
            else params.mismatch_rate / 3.0
        )

    @lru_cache(maxsize=None)
    def visit(state: str, i: int, j: int) -> float:
        if i == n:
            # Termination is free after all read bases have been emitted.
            # Deletion chains after the final emission are separate valid
            # paths represented in the forward row and must also be counted.
            total = 1.0 if j >= 1 else 0.0
            if j < m and state in {"M", "D"}:
                _, _, md, _, gap_extend = transition(j)
                previous = md if state == "M" else gap_extend
                total += previous * visit("D", i, j + 1)
            return total
        total = 0.0
        if j < m:
            mm, _, _, gap_to_m, _ = transition(j)
            previous = mm if state == "M" else gap_to_m
            total += (
                previous
                * match_emission(read[i], template[j])
                * visit("M", i + 1, j + 1)
            )
        if j >= 1:
            _, mi, _, _, gap_extend = transition(j - 1)
            previous = mi if state == "M" else gap_extend if state == "I" else 0
            total += (
                previous
                * params.insertion_background
                * visit("I", i + 1, j)
            )
        if i > 0 and j < m:
            _, _, md, _, gap_extend = transition(j)
            previous = md if state == "M" else gap_extend if state == "D" else 0
            total += previous * visit("D", i, j + 1)
        return total

    return sum(
        visit("M", 0, start) / (m + 1.0) for start in range(m + 1)
    )


def test_phmm_forward_matches_bruteforce() -> None:
    params = PHMMParameters(
        mismatch_rate=0.02,
        flank_indel_open=0.04,
        repeat_indel_open=0.08,
        flank_indel_extend=0.20,
        repeat_indel_extend=0.40,
    )
    model = MotifAwarePHMM(params)
    candidate = allele("ACG", 1, 2)
    expected = brute_force_probability("AG", candidate, params)
    observed = math.exp(model._forward_log_likelihood_python("AG", candidate))
    assert math.isclose(observed, expected, rel_tol=1e-12, abs_tol=1e-15)
