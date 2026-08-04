from __future__ import annotations

import math
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator

from .index import SeedIndex
from .io import iter_fastq
from .models import MappingAssignment, SeedHit
from .syncmer import iter_syncmers
from .utils import open_text


@dataclass
class CandidateChain:
    pointer_id: str
    hits: list[SeedHit]
    score: float
    path_start: int
    path_end: int
    read_start: int
    read_end: int
    left_hits: int
    right_hits: int
    support_fraction: float
    dispersion: float


def _best_chain_for_locus(
    pointer_id: str,
    hits: list[SeedHit],
    flank_seed_count: int,
    k: int,
    delta: int,
    tau0: float,
    min_seeds: int,
) -> CandidateChain | None:
    bins: dict[int, list[SeedHit]] = defaultdict(list)
    for hit in hits:
        diagonal = hit.path_pos - hit.read_pos
        bins[round(diagonal / delta)].append(hit)
    best_hits = max(bins.values(), key=lambda values: (len(values), max(h.path_pos for h in values) - min(h.path_pos for h in values)))
    unique = {(h.seed, h.read_pos, h.path_pos, h.side): h for h in best_hits}
    chain_hits = sorted(unique.values(), key=lambda h: (h.path_pos, h.read_pos))
    if len(chain_hits) < min_seeds:
        return None
    diagonals = [h.path_pos - h.read_pos for h in chain_hits]
    mean_d = sum(diagonals) / len(diagonals)
    dispersion = math.sqrt(sum((x - mean_d) ** 2 for x in diagonals) / len(diagonals))
    support_fraction = len(chain_hits) / max(1, flank_seed_count)
    threshold = max(tau0, 2.0 / max(1, flank_seed_count))
    if support_fraction < threshold:
        return None
    left_hits = sum(h.side == "L" for h in chain_hits)
    right_hits = sum(h.side == "R" for h in chain_hits)
    span = max(h.path_pos for h in chain_hits) - min(h.path_pos for h in chain_hits) + k
    score = len(chain_hits) + 0.01 * span - 0.02 * dispersion
    return CandidateChain(
        pointer_id=pointer_id,
        hits=chain_hits,
        score=score,
        path_start=min(h.path_pos for h in chain_hits),
        path_end=max(h.path_pos for h in chain_hits) + k,
        read_start=min(h.read_pos for h in chain_hits),
        read_end=max(h.read_pos for h in chain_hits) + k,
        left_hits=left_hits,
        right_hits=right_hits,
        support_fraction=support_fraction,
        dispersion=dispersion,
    )


def _coarse_seed_count_for_locus(
    pointer_id: str,
    hits: list[SeedHit],
    flank_seed_count: int,
    k: int,
    tau0: float,
    min_seeds: int,
) -> CandidateChain | None:
    """No-chaining control based only on unique seed-hit count by pointer."""
    unique = {(h.seed, h.read_pos, h.path_pos, h.side): h for h in hits}
    coarse_hits = list(unique.values())
    if len(coarse_hits) < min_seeds:
        return None
    support_fraction = len(coarse_hits) / max(1, flank_seed_count)
    threshold = max(tau0, 2.0 / max(1, flank_seed_count))
    if support_fraction < threshold:
        return None
    return CandidateChain(
        pointer_id=pointer_id,
        hits=coarse_hits,
        score=float(len(coarse_hits)),
        path_start=min(h.path_pos for h in coarse_hits),
        path_end=max(h.path_pos for h in coarse_hits) + k,
        read_start=min(h.read_pos for h in coarse_hits),
        read_end=max(h.read_pos for h in coarse_hits) + k,
        left_hits=sum(h.side == "L" for h in coarse_hits),
        right_hits=sum(h.side == "R" for h in coarse_hits),
        support_fraction=support_fraction,
        dispersion=0.0,
    )


def map_sequence(
    seq: str,
    index: SeedIndex,
    delta: int = 10,
    tau0: float = 0.20,
    min_seeds: int = 2,
    max_seed_occurrences: int = 200,
    use_chaining: bool = True,
) -> list[CandidateChain]:
    k = int(index.metadata["k"])
    s = int(index.metadata["s"])
    t = int(index.metadata["t"])
    by_locus: dict[str, list[SeedHit]] = defaultdict(list)
    for syncmer in iter_syncmers(seq, k=k, s=s, t=t):
        for row in index.lookup(syncmer.seed, max_occurrences=max_seed_occurrences):
            by_locus[row["pointer_id"]].append(
                SeedHit(
                    pointer_id=row["pointer_id"],
                    side=row["side"],
                    path_pos=int(row["path_pos"]),
                    read_pos=syncmer.position,
                    orientation=syncmer.orientation * int(row["orientation"]),
                    seed=syncmer.seed,
                )
            )
    chains: list[CandidateChain] = []
    for pointer_id, hits in by_locus.items():
        locus_row = index.locus(pointer_id)
        if use_chaining:
            chain = _best_chain_for_locus(
                pointer_id,
                hits,
                int(locus_row["flank_seed_count"]),
                k,
                delta,
                tau0,
                min_seeds,
            )
        else:
            chain = _coarse_seed_count_for_locus(
                pointer_id,
                hits,
                int(locus_row["flank_seed_count"]),
                k,
                tau0,
                min_seeds,
            )
        if chain:
            chains.append(chain)
    return sorted(chains, key=lambda c: (c.score, len(c.hits), -c.dispersion), reverse=True)


def _combine_pair(c1: list[CandidateChain], c2: list[CandidateChain]) -> list[CandidateChain]:
    by1 = {c.pointer_id: c for c in c1}
    by2 = {c.pointer_id: c for c in c2}
    combined: list[CandidateChain] = []
    for pid in set(by1) | set(by2):
        a, b = by1.get(pid), by2.get(pid)
        source = a or b
        assert source is not None
        if a and b:
            combined.append(
                CandidateChain(
                    pointer_id=pid,
                    hits=a.hits + b.hits,
                    score=a.score + b.score + 0.5,
                    path_start=min(a.path_start, b.path_start),
                    path_end=max(a.path_end, b.path_end),
                    read_start=min(a.read_start, b.read_start),
                    read_end=max(a.read_end, b.read_end),
                    left_hits=a.left_hits + b.left_hits,
                    right_hits=a.right_hits + b.right_hits,
                    support_fraction=max(a.support_fraction, b.support_fraction),
                    dispersion=(a.dispersion + b.dispersion) / 2,
                )
            )
        else:
            combined.append(source)
    return sorted(combined, key=lambda c: c.score, reverse=True)


def _mapq(best: float, second: float | None) -> int:
    if second is None:
        return 60
    margin = max(0.0, best - second)
    return min(60, int(round(10 * margin)))


def _format_gaf(name: str, seq_len: int, chain: CandidateChain, path_length: int, mapq: int, mate: int) -> str:
    matches = len(chain.hits)
    block = max(1, chain.read_end - chain.read_start)
    tags = [
        f"pi:Z:{chain.pointer_id}",
        f"as:f:{chain.score:.4f}",
        f"fc:i:{len(chain.hits)}",
        f"lh:i:{chain.left_hits}",
        f"rh:i:{chain.right_hits}",
        f"sf:f:{chain.support_fraction:.6f}",
        f"ds:f:{chain.dispersion:.6f}",
        f"mt:i:{mate}",
    ]
    return "\t".join(
        [
            name,
            str(seq_len),
            str(chain.read_start),
            str(min(seq_len, chain.read_end)),
            "+",
            chain.pointer_id,
            str(path_length),
            str(chain.path_start),
            str(min(path_length, chain.path_end)),
            str(matches),
            str(block),
            str(mapq),
            *tags,
        ]
    )


def map_fastq(
    index_path: str | Path,
    reads1: str | Path,
    out_path: str | Path,
    reads2: str | Path | None = None,
    delta: int = 10,
    tau0: float = 0.20,
    min_seeds: int = 2,
    min_score_margin: float = 0.25,
    max_seed_occurrences: int = 200,
    use_chaining: bool = True,
) -> dict[str, int]:
    index = SeedIndex(index_path)
    stats = {"pairs": 0, "mapped": 0, "ambiguous": 0, "unmapped": 0}
    fq2_iter = iter_fastq(reads2) if reads2 else None
    with open_text(out_path, "wt") as out:
        for rec1 in iter_fastq(reads1):
            rec2 = next(fq2_iter) if fq2_iter else None
            if rec2 and rec1[0] != rec2[0]:
                raise ValueError(f"Paired FASTQ names differ: {rec1[0]} vs {rec2[0]}")
            name, seq1, _ = rec1
            c1 = map_sequence(
                seq1, index, delta, tau0, min_seeds, max_seed_occurrences, use_chaining
            )
            c2 = (
                map_sequence(
                    rec2[1],
                    index,
                    delta,
                    tau0,
                    min_seeds,
                    max_seed_occurrences,
                    use_chaining,
                )
                if rec2
                else []
            )
            candidates = _combine_pair(c1, c2) if rec2 else c1
            stats["pairs"] += 1
            if not candidates:
                stats["unmapped"] += 1
                continue
            best = candidates[0]
            second_score = candidates[1].score if len(candidates) > 1 else None
            if second_score is not None and best.score - second_score < min_score_margin:
                stats["ambiguous"] += 1
                continue
            row = index.locus(best.pointer_id)
            mq = _mapq(best.score, second_score)
            out.write(_format_gaf(name, len(seq1), best, int(row["path_length"]), mq, 1 if rec2 else 0) + "\n")
            if rec2:
                out.write(_format_gaf(name, len(rec2[1]), best, int(row["path_length"]), mq, 2) + "\n")
            stats["mapped"] += 1
    index.close()
    return stats


def map_bam(
    index_path: str | Path,
    bam_path: str | Path,
    out_path: str | Path,
    delta: int = 10,
    tau0: float = 0.20,
    min_seeds: int = 2,
    min_score_margin: float = 0.25,
    max_seed_occurrences: int = 200,
    use_chaining: bool = True,
) -> dict[str, int]:
    from .io import iter_bam_reads

    index = SeedIndex(index_path)
    stats = {"reads": 0, "mapped": 0, "ambiguous": 0, "unmapped": 0}
    with open_text(out_path, "wt") as out:
        for name, seq, mate in iter_bam_reads(bam_path):
            candidates = map_sequence(
                seq,
                index,
                delta,
                tau0,
                min_seeds,
                max_seed_occurrences,
                use_chaining,
            )
            stats["reads"] += 1
            if not candidates:
                stats["unmapped"] += 1
                continue
            best = candidates[0]
            second_score = candidates[1].score if len(candidates) > 1 else None
            if second_score is not None and best.score - second_score < min_score_margin:
                stats["ambiguous"] += 1
                continue
            row = index.locus(best.pointer_id)
            mq = _mapq(best.score, second_score)
            out.write(_format_gaf(name, len(seq), best, int(row["path_length"]), mq, mate) + "\n")
            stats["mapped"] += 1
    index.close()
    return stats
