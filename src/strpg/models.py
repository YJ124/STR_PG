from __future__ import annotations

from dataclasses import dataclass, field, asdict
from typing import Any


@dataclass
class AlleleRecord:
    allele_id: str
    repeat_count: int
    sequence: str
    repeat_sequence: str
    motif: str
    repeat_start: int
    repeat_end: int
    frequencies: dict[str, float] = field(default_factory=dict)
    counts: dict[str, float] = field(default_factory=dict)
    status: str = "registered"
    source: str = "catalog"
    metadata: dict[str, Any] = field(default_factory=dict)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "AlleleRecord":
        return cls(**data)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class ProvisionalAllele:
    repeat_count: int
    sequence: str
    support: int
    flank_support: int
    log_likelihood_gain: float
    gq: int
    status: str
    reason: str
    source_sample: str | None = None
    metadata: dict[str, Any] = field(default_factory=dict)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "ProvisionalAllele":
        return cls(**data)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class LocusRecord:
    pointer_id: str
    build: str
    chrom: str
    start: int
    end: int
    motif: str
    left_flank: str
    right_flank: str
    anchor_length: int
    alleles: list[AlleleRecord] = field(default_factory=list)
    provisional: list[ProvisionalAllele] = field(default_factory=list)
    metadata: dict[str, Any] = field(default_factory=dict)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "LocusRecord":
        payload = dict(data)
        payload["alleles"] = [AlleleRecord.from_dict(x) for x in payload.get("alleles", [])]
        payload["provisional"] = [ProvisionalAllele.from_dict(x) for x in payload.get("provisional", [])]
        return cls(**payload)

    def to_dict(self) -> dict[str, Any]:
        return {
            **asdict(self),
            "alleles": [a.to_dict() for a in self.alleles],
            "provisional": [p.to_dict() for p in self.provisional],
        }

    def allele_by_count(self, repeat_count: int) -> AlleleRecord | None:
        for allele in self.alleles:
            if allele.repeat_count == repeat_count:
                return allele
        return None


@dataclass
class SeedHit:
    pointer_id: str
    side: str
    path_pos: int
    read_pos: int
    orientation: int
    seed: str


@dataclass
class MappingAssignment:
    query_name: str
    pointer_id: str
    score: float
    mapq: int
    read_start: int
    read_end: int
    path_start: int
    path_end: int
    path_length: int
    strand: str = "+"
    mate: int = 0
    tags: dict[str, str | int | float] = field(default_factory=dict)


@dataclass
class GenotypeCall:
    pointer_id: str
    chrom: str
    start: int
    end: int
    motif: str
    depth: int
    allele1: int
    allele2: int
    gt: str
    gq: int
    posterior: float
    log_likelihood: float
    log_prior: float
    log_score: float
    log_mix: float | None
    d_paper: float | None
    d_likelihood: float | None
    mixture_weights: dict[int, float] | None
    heterogeneity_threshold: float | None
    heterogeneity_flag: str
