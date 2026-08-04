from __future__ import annotations

from collections import Counter

from .models import LocusRecord
from .utils import reverse_complement


def has_anchored_repeat_evidence(
    read: str,
    locus: LocusRecord,
    anchor_k: int = 12,
    min_repeat_units: int = 2,
    max_repeat_mismatch_fraction: float = 0.10,
) -> bool:
    """Return whether a read connects a unique flank anchor to STR sequence.

    Pure-flank reads are useful for locus localization but are not informative
    about allele length.  Requiring an observed flank/repeat boundary prevents
    such reads from contributing allele-dependent pHMM path multiplicity.
    """
    if len(locus.left_flank) < anchor_k or len(locus.right_flank) < anchor_k:
        return False
    motif = locus.motif.upper()
    minimum = len(motif) * min_repeat_units
    left_anchor = locus.left_flank[-anchor_k:].upper()
    right_anchor = locus.right_flank[:anchor_k].upper()
    for oriented in (read.upper(), reverse_complement(read.upper())):
        left_pos = oriented.find(left_anchor)
        if left_pos >= 0:
            repeat_start = left_pos + anchor_k
            observed = oriented[repeat_start:repeat_start + minimum]
            if len(observed) == minimum and _motif_mismatch_fraction(observed, motif) <= max_repeat_mismatch_fraction:
                return True
        right_pos = oriented.find(right_anchor)
        if right_pos >= minimum:
            observed = oriented[right_pos - minimum:right_pos]
            if _motif_mismatch_fraction(observed, motif) <= max_repeat_mismatch_fraction:
                return True
    return False


def _motif_mismatch_fraction(sequence: str, motif: str) -> float:
    if not sequence:
        return 1.0
    expected = (motif * ((len(sequence) + len(motif) - 1) // len(motif)))[: len(sequence)]
    return sum(a != b for a, b in zip(sequence.upper(), expected.upper())) / len(sequence)


def infer_spanning_repeat_count(
    read: str,
    locus: LocusRecord,
    anchor_k: int = 12,
    max_repeat_mismatch_fraction: float = 0.10,
) -> int | None:
    """Infer a repeat count only from reads containing unique anchors on both sides.

    This deliberately conservative rule implements the manuscript requirement that
    novel alleles be supported by flanking evidence rather than repeat sequence alone.
    """
    if len(locus.left_flank) < anchor_k or len(locus.right_flank) < anchor_k:
        return None
    left_anchor = locus.left_flank[-anchor_k:].upper()
    right_anchor = locus.right_flank[:anchor_k].upper()
    for oriented in (read.upper(), reverse_complement(read.upper())):
        left_pos = oriented.find(left_anchor)
        if left_pos < 0:
            continue
        repeat_start = left_pos + anchor_k
        right_pos = oriented.find(right_anchor, repeat_start)
        if right_pos < 0:
            continue
        repeat_sequence = oriented[repeat_start:right_pos]
        if len(repeat_sequence) % len(locus.motif) != 0:
            continue
        if _motif_mismatch_fraction(repeat_sequence, locus.motif) > max_repeat_mismatch_fraction:
            continue
        return len(repeat_sequence) // len(locus.motif)
    return None


def candidate_repeat_counts(
    reads: list[str],
    locus: LocusRecord,
    min_support: int = 3,
    anchor_k: int = 12,
    max_repeat_mismatch_fraction: float = 0.10,
) -> Counter[int]:
    counts: Counter[int] = Counter()
    registered = {a.repeat_count for a in locus.alleles}
    for read in reads:
        repeat_count = infer_spanning_repeat_count(read, locus, anchor_k, max_repeat_mismatch_fraction)
        if repeat_count is not None and repeat_count not in registered:
            counts[repeat_count] += 1
    return Counter({k: v for k, v in counts.items() if v >= min_support})
