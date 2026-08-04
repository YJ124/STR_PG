from strpg.models import AlleleRecord
from strpg.phmm import MotifAwarePHMM


def make(count):
    left="ACGTTGCAACGT"
    right="TGCACGTTTGCA"
    seq=left+"CAG"*count+right
    return AlleleRecord(f"L{count}",count,seq,"CAG"*count,"CAG",len(left),len(left)+3*count,{"GLOBAL":1.0},{"GLOBAL":1000})


def test_matching_allele_scores_higher():
    hmm=MotifAwarePHMM()
    a8,a12=make(8),make(12)
    read=a12.sequence
    assert hmm.forward_log_likelihood(read,a12) > hmm.forward_log_likelihood(read,a8)
