from strpg.models import AlleleRecord, LocusRecord
from strpg.novel import has_anchored_repeat_evidence
from strpg.phmm import PHMMParameters


def make_locus() -> LocusRecord:
    left = "AACCGGTTAACCGGTTAACC"
    right = "TTGGCCAATTGGCCAATTGG"
    motif = "CAG"
    allele = AlleleRecord(
        allele_id="P1:L10",
        repeat_count=10,
        sequence=left + motif * 10 + right,
        repeat_sequence=motif * 10,
        motif=motif,
        repeat_start=len(left),
        repeat_end=len(left) + 30,
        frequencies={"GLOBAL": 1.0},
        counts={"GLOBAL": 1000.0},
    )
    return LocusRecord(
        pointer_id="P1",
        build="GRCh38",
        chrom="chr1",
        start=100,
        end=129,
        motif=motif,
        left_flank=left,
        right_flank=right,
        anchor_length=len(left),
        alleles=[allele],
    )


def test_flank_repeat_boundary_is_informative_but_pure_flank_is_not():
    locus = make_locus()
    boundary_read = locus.left_flank[-12:] + locus.motif * 4
    pure_flank = locus.left_flank
    assert has_anchored_repeat_evidence(boundary_read, locus)
    assert not has_anchored_repeat_evidence(pure_flank, locus)


def test_calibrated_repeat_transition_defaults():
    params = PHMMParameters()
    assert params.repeat_indel_open == 0.005
    assert params.repeat_indel_extend == 0.30
