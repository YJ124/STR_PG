from strpg.models import AlleleRecord
from strpg.priors import balding_nichols_prior, genotype_priors, length_similarity


def allele(count, freq):
    return AlleleRecord(str(count), count, "A", "A", "A", 0, 1, {"GLOBAL": freq}, {"GLOBAL": freq*1000})


def test_bn_and_smoothing():
    assert balding_nichols_prior(0.5, 0.5, True, 0.1) == 0.275
    assert length_similarity(10, 10, 2) == 1.0
    assert length_similarity(10, 20, 2) < 1e-4
    priors = genotype_priors([allele(10, 0.5), allele(12, 0.5)], {"GLOBAL": 1.0}, theta=0.01, sigma=10)
    assert abs(sum(priors.values()) - 1) < 1e-12
