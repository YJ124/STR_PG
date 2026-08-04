import numpy as np

from strpg.genotyper import call_locus
from strpg.phmm import MotifAwarePHMM


def test_mixture_does_not_change_gt_gq(read_selection_fixture):
    locus = read_selection_fixture["locus"]
    reads = ["A", "C"]
    matrix = np.array([[-1.0, -5.0, -8.0], [-1.0, -5.0, -8.0]])
    kwargs = dict(
        locus=locus, reads=reads, phmm=MotifAwarePHMM(),
        population_mix={"GLOBAL": 1.0}, theta=0.01, sigma=10.0,
        precomputed_log_likelihoods=matrix,
    )
    with_mix, _, _ = call_locus(**kwargs, use_mixture_model=True)
    without_mix, _, _ = call_locus(**kwargs, use_mixture_model=False)
    assert (with_mix.gt, with_mix.gq) == (without_mix.gt, without_mix.gq)
