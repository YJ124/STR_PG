import math
import random

from conftest import make_locus
from strpg.phmm import MotifAwarePHMM


def test_jit_phmm_matches_python_recurrence_to_numerical_precision():
    rng = random.Random(7242026)
    phmm = MotifAwarePHMM()
    locus = make_locus()
    for allele in locus.alleles:
        reads = [
            allele.sequence[5:55],
            allele.sequence[-55:-5],
            "".join(rng.choice("ACGT") for _ in range(70)),
        ]
        for read in reads:
            expected = phmm._forward_log_likelihood_python(read, allele)
            observed = phmm.forward_log_likelihood(read, allele)
            assert math.isclose(observed, expected, rel_tol=0.0, abs_tol=1e-12)
