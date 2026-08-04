from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Iterable

from .models import AlleleRecord, LocusRecord, ProvisionalAllele
from .utils import open_text


class AlleleRegistry:
    def __init__(self, loci: dict[str, LocusRecord] | None = None):
        self.loci = loci or {}

    @classmethod
    def load(cls, path: str | Path) -> "AlleleRegistry":
        loci: dict[str, LocusRecord] = {}
        with open_text(path, "rt") as handle:
            for raw in handle:
                if not raw.strip():
                    continue
                locus = LocusRecord.from_dict(json.loads(raw))
                if locus.pointer_id in loci:
                    raise ValueError(f"Duplicate pointer_id in registry: {locus.pointer_id}")
                loci[locus.pointer_id] = locus
        return cls(loci)

    def save(self, path: str | Path) -> None:
        tmp = Path(str(path) + ".tmp")
        with open_text(tmp, "wt") as handle:
            for pointer_id in sorted(self.loci):
                handle.write(json.dumps(self.loci[pointer_id].to_dict(), ensure_ascii=False, sort_keys=True) + "\n")
        tmp.replace(path)

    def validate(self) -> list[str]:
        errors: list[str] = []
        for pointer_id, locus in self.loci.items():
            if not locus.alleles:
                errors.append(f"{pointer_id}: no registered alleles")
            seen: set[int] = set()
            for allele in locus.alleles:
                if allele.repeat_count in seen:
                    errors.append(f"{pointer_id}: duplicate repeat count {allele.repeat_count}")
                seen.add(allele.repeat_count)
                expected = locus.left_flank + locus.motif * allele.repeat_count + locus.right_flank
                if allele.sequence != expected:
                    errors.append(f"{pointer_id}/{allele.allele_id}: sequence does not equal flank+motif*L+flank")
                if allele.repeat_start != len(locus.left_flank):
                    errors.append(f"{pointer_id}/{allele.allele_id}: invalid repeat_start")
                if allele.repeat_end != len(locus.left_flank) + len(locus.motif) * allele.repeat_count:
                    errors.append(f"{pointer_id}/{allele.allele_id}: invalid repeat_end")
        return errors

    def populations(self) -> list[str]:
        pops: set[str] = set()
        for locus in self.loci.values():
            for allele in locus.alleles:
                pops.update(allele.frequencies)
                pops.update(allele.counts)
        return sorted(pops)

    @staticmethod
    def recompute_frequencies(locus: LocusRecord, populations: Iterable[str] | None = None) -> None:
        pops = set(populations or [])
        for allele in locus.alleles:
            pops.update(allele.counts)
            pops.update(allele.frequencies)
        for pop in pops:
            total = sum(max(0.0, a.counts.get(pop, 0.0)) for a in locus.alleles)
            if total <= 0:
                continue
            for allele in locus.alleles:
                allele.frequencies[pop] = max(0.0, allele.counts.get(pop, 0.0)) / total

    @staticmethod
    def initialize_counts_from_frequencies(locus: LocusRecord, base_total: float = 1000.0) -> None:
        pops: set[str] = set()
        for allele in locus.alleles:
            pops.update(allele.frequencies)
        for pop in pops:
            for allele in locus.alleles:
                allele.counts.setdefault(pop, allele.frequencies.get(pop, 0.0) * base_total)

    def add_provisional(self, pointer_id: str, provisional: ProvisionalAllele) -> None:
        locus = self.loci[pointer_id]
        for idx, old in enumerate(locus.provisional):
            if old.repeat_count == provisional.repeat_count:
                if provisional.support >= old.support:
                    locus.provisional[idx] = provisional
                return
        locus.provisional.append(provisional)

    def accept_novel(
        self,
        pointer_id: str,
        repeat_count: int,
        population: str,
        alpha: float = 0.5,
        source_sample: str | None = None,
    ) -> AlleleRecord:
        locus = self.loci[pointer_id]
        existing = locus.allele_by_count(repeat_count)
        if existing:
            return existing
        populations = set(self.populations()) | {population}
        sequence = locus.left_flank + locus.motif * repeat_count + locus.right_flank
        allele = AlleleRecord(
            allele_id=f"{pointer_id}:L{repeat_count}",
            repeat_count=repeat_count,
            sequence=sequence,
            repeat_sequence=locus.motif * repeat_count,
            motif=locus.motif,
            repeat_start=len(locus.left_flank),
            repeat_end=len(locus.left_flank) + len(locus.motif) * repeat_count,
            frequencies={},
            counts={pop: alpha + (1.0 if pop == population else 0.0) for pop in populations},
            status="registered",
            source="novel_discovery",
            metadata={"source_sample": source_sample},
        )
        for old in locus.alleles:
            for pop in populations:
                old.counts.setdefault(pop, alpha)
        locus.alleles.append(allele)
        locus.alleles.sort(key=lambda a: a.repeat_count)
        self.recompute_frequencies(locus, populations)
        locus.provisional = [p for p in locus.provisional if p.repeat_count != repeat_count]
        return allele
