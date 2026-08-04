from __future__ import annotations

from dataclasses import dataclass
from typing import Iterator

from .utils import canonical_kmer, stable_hash


@dataclass(frozen=True)
class Syncmer:
    seed: str
    position: int
    orientation: int


def iter_syncmers(seq: str, k: int = 15, s: int = 5, t: int = 2) -> Iterator[Syncmer]:
    """Yield deterministic open syncmers.

    A k-mer is selected when the minimum-hash s-mer within it begins at offset t.
    The yielded seed is canonicalized so both strands share an index key.
    """
    seq = seq.upper()
    if not (1 <= s <= k):
        raise ValueError("Require 1 <= s <= k")
    if not (0 <= t <= k - s):
        raise ValueError("Require 0 <= t <= k-s")
    # Each s-mer belongs to as many as k-s+1 overlapping k-mers.  Hash every
    # s-mer once, then reuse the exact integer hashes in each window.  This is
    # mathematically identical to the direct definition above and leaves the
    # first-minimum tie rule unchanged.
    s_hashes = [
        stable_hash(seq[i : i + s])
        if all(base in "ACGT" for base in seq[i : i + s])
        else None
        for i in range(max(0, len(seq) - s + 1))
    ]
    window_width = k - s + 1
    for i in range(0, len(seq) - k + 1):
        kmer = seq[i : i + k]
        if any(base not in "ACGT" for base in kmer):
            continue
        window = s_hashes[i : i + window_width]
        if min(range(window_width), key=window.__getitem__) != t:
            continue
        canonical, orientation = canonical_kmer(kmer)
        yield Syncmer(canonical, i, orientation)
