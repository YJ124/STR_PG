from strpg.models import AlleleRecord
from strpg.phmm import MotifAwarePHMM


def allele(count):
    left, right, motif = "ACGT" * 10, "TGCA" * 10, "CAG"
    repeat = motif * count
    return AlleleRecord(
        f"L{count}", count, left + repeat + right, repeat, motif,
        len(left), len(left) + len(repeat),
    )


def test_phmm_exact_match_ranking():
    candidates = [allele(8), allele(12), allele(16)]
    read = candidates[1].sequence
    scores = [
        MotifAwarePHMM().orientation_marginal_log_likelihood(read, candidate)
        for candidate in candidates
    ]
    assert scores.index(max(scores)) == 1
