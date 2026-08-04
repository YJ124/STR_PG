import numpy as np
from strpg.mixture import fit_allele_mixture


def test_em_recovers_three_components():
    rows=[]
    for k in range(3):
        for _ in range(20):
            row=np.full(3,-20.0)
            row[k]=0.0
            rows.append(row)
    fit=fit_allele_mixture(np.asarray(rows))
    assert fit.converged
    assert np.allclose(fit.weights,[1/3]*3,atol=0.02)
