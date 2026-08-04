"""Deterministic, likelihood-backend-independent informative-read selection."""

from __future__ import annotations

import csv
import hashlib
import re
from dataclasses import asdict, dataclass, replace
from pathlib import Path

from .defaults import READ_CACHE_SCHEMA_VERSION
from .io import iter_fastq, read_gaf
from .models import LocusRecord
from .novel import has_anchored_repeat_evidence
from .registry import AlleleRegistry
from .utils import reverse_complement


def _sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("ascii")).hexdigest()


@dataclass(frozen=True)
class ReadCacheRecord:
    schema_version: str
    sample_id: str
    locus_id: str
    read_id: str
    mate: int
    sequence: str
    base_quality: str
    mapq: int
    selected: bool
    selection_reason: str
    left_anchor: bool
    right_anchor: bool
    repeat_evidence: bool
    sequence_sha256: str
    source_fastq_sha256: str
    source_gaf_sha256: str


READ_CACHE_FIELDS = list(ReadCacheRecord.__dataclass_fields__)


def _evidence(sequence: str, locus: LocusRecord, anchor_k: int) -> tuple[bool, bool, bool]:
    orientations = (sequence.upper(), reverse_complement(sequence.upper()))
    left = any(
        locus.left_flank[-anchor_k:].upper() in oriented
        for oriented in orientations
    )
    right = any(
        locus.right_flank[:anchor_k].upper() in oriented
        for oriented in orientations
    )
    repeat = any(
        locus.motif.upper() * 2 in oriented for oriented in orientations
    )
    anchored = has_anchored_repeat_evidence(
        sequence, locus, anchor_k=anchor_k
    )
    return left, right, repeat or anchored


def select_informative_read_records(
    *,
    sample_id: str,
    registry: AlleleRegistry,
    gaf_path: str | Path,
    reads1: str | Path,
    reads2: str | Path | None = None,
    min_mapq: int = 0,
    anchor_k: int = 12,
) -> list[ReadCacheRecord]:
    """Return deterministic records without constructing a likelihood backend."""

    assignments: dict[tuple[str, int], tuple[str, int]] = {}
    for alignment in read_gaf(gaf_path):
        if alignment.mapq >= min_mapq:
            assignments[(alignment.query_name, alignment.mate)] = (
                alignment.pointer_id,
                alignment.mapq,
            )
    gaf_hash = _sha256_file(gaf_path)
    fastq_hashes = {1: _sha256_file(reads1)}
    if reads2 is not None:
        fastq_hashes[2] = _sha256_file(reads2)

    records: list[ReadCacheRecord] = []
    iterators = [(1, iter_fastq(reads1))]
    if reads2 is not None:
        iterators.append((2, iter_fastq(reads2)))
    for mate, iterator in iterators:
        for read_id, sequence, quality in iterator:
            assignment = assignments.get((read_id, mate)) or assignments.get(
                (read_id, 0)
            )
            if assignment is None:
                continue
            locus_id, mapq = assignment
            locus = registry.loci.get(locus_id)
            if locus is None:
                continue
            left, right, repeat = _evidence(sequence, locus, anchor_k)
            selected = has_anchored_repeat_evidence(
                sequence, locus, anchor_k=anchor_k
            )
            reason = (
                "FLANK_REPEAT_BOUNDARY"
                if selected
                else "MAPPED_WITHOUT_BOUNDARY_EVIDENCE"
            )
            records.append(
                ReadCacheRecord(
                    schema_version=READ_CACHE_SCHEMA_VERSION,
                    sample_id=sample_id,
                    locus_id=locus_id,
                    read_id=read_id,
                    mate=mate,
                    sequence=sequence,
                    base_quality=quality,
                    mapq=mapq,
                    selected=selected,
                    selection_reason=reason,
                    left_anchor=left,
                    right_anchor=right,
                    repeat_evidence=repeat,
                    sequence_sha256=_sha256_text(sequence),
                    source_fastq_sha256=fastq_hashes[mate],
                    source_gaf_sha256=gaf_hash,
                )
            )
    by_pair: dict[tuple[str, str], list[int]] = {}
    for index, record in enumerate(records):
        by_pair.setdefault((record.locus_id, record.read_id), []).append(index)
    for indices in by_pair.values():
        if len(indices) > 1 and any(records[index].selected for index in indices):
            for index in indices:
                records[index] = replace(
                    records[index],
                    selected=True,
                    selection_reason="PAIRED_WITH_BOUNDARY_EVIDENCE",
                )
    return sorted(records, key=lambda r: (r.locus_id, r.read_id, r.mate))


def write_read_cache(records: list[ReadCacheRecord], cache_root: str | Path) -> list[Path]:
    """Write byte-stable per-sample/per-locus TSV files."""

    root = Path(cache_root)
    grouped: dict[tuple[str, str], list[ReadCacheRecord]] = {}
    for record in records:
        grouped.setdefault((record.sample_id, record.locus_id), []).append(record)
    paths: list[Path] = []
    for (sample_id, locus_id), group in sorted(grouped.items()):
        safe_sample = re.sub(r"[^A-Za-z0-9._-]+", "_", sample_id).strip("._")
        if not safe_sample:
            raise ValueError("sample_id does not contain a filesystem-safe character")
        safe_locus = re.sub(r"[^A-Za-z0-9._-]+", "_", locus_id).strip("._")
        if not safe_locus:
            safe_locus = "locus"
        if safe_locus != locus_id:
            safe_locus += "_" + hashlib.sha256(
                locus_id.encode("utf-8")
            ).hexdigest()[:8]
        path = root / safe_sample / f"{safe_locus}.tsv"
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=READ_CACHE_FIELDS,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            for record in group:
                writer.writerow(asdict(record))
        paths.append(path)
    return paths


def selected_reads_by_locus(
    records: list[ReadCacheRecord],
) -> tuple[dict[str, list[str]], dict[str, list[str]], dict[str, list[str]]]:
    """Return mapped sequences, selected sequences, and selected read IDs."""

    mapped: dict[str, list[str]] = {}
    selected: dict[str, list[str]] = {}
    selected_ids: dict[str, list[str]] = {}
    for record in records:
        mapped.setdefault(record.locus_id, []).append(record.sequence)
        if record.selected:
            selected.setdefault(record.locus_id, []).append(record.sequence)
            selected_ids.setdefault(record.locus_id, []).append(
                f"{record.read_id}/{record.mate}"
            )
    return mapped, selected, selected_ids
