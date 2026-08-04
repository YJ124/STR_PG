from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Iterator

from .models import MappingAssignment
from .utils import open_text


def read_fasta(path: str | Path) -> dict[str, str]:
    sequences: dict[str, list[str]] = {}
    name: str | None = None
    with open_text(path, "rt") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].split()[0]
                if name in sequences:
                    raise ValueError(f"Duplicate FASTA record: {name}")
                sequences[name] = []
            elif name is None:
                raise ValueError("FASTA sequence encountered before header")
            else:
                sequences[name].append(line.upper())
    return {k: "".join(v) for k, v in sequences.items()}


def iter_fastq(path: str | Path) -> Iterator[tuple[str, str, str]]:
    with open_text(path, "rt") as handle:
        while True:
            h = handle.readline()
            if not h:
                break
            seq = handle.readline().strip()
            plus = handle.readline()
            qual = handle.readline().strip()
            if not plus or not qual:
                raise ValueError(f"Truncated FASTQ: {path}")
            name = h.strip()[1:].split()[0]
            if name.endswith("/1") or name.endswith("/2"):
                name = name[:-2]
            yield name, seq.upper(), qual


def write_fastq(records: Iterator[tuple[str, str, str]], path: str | Path) -> None:
    with open_text(path, "wt") as handle:
        for name, seq, qual in records:
            handle.write(f"@{name}\n{seq}\n+\n{qual}\n")


def read_catalog(path: str | Path) -> list[dict[str, str]]:
    with open_text(path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"build", "chrom", "start", "end", "motif", "candidate_Ls", "pointer_id"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"STR catalog missing columns: {sorted(missing)}")
        return [dict(row) for row in reader]


def read_gaf(path: str | Path) -> Iterator[MappingAssignment]:
    with open_text(path, "rt") as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 12:
                raise ValueError(f"Invalid GAF line with {len(fields)} fields")
            tags: dict[str, str | int | float] = {}
            for tag in fields[12:]:
                parts = tag.split(":", 2)
                if len(parts) != 3:
                    continue
                key, typ, value = parts
                if typ == "i":
                    tags[key] = int(value)
                elif typ == "f":
                    tags[key] = float(value)
                else:
                    tags[key] = value
            yield MappingAssignment(
                query_name=fields[0],
                pointer_id=str(tags.get("pi", fields[5])),
                score=float(tags.get("as", fields[9])),
                mapq=int(fields[11]),
                read_start=int(fields[2]),
                read_end=int(fields[3]),
                path_start=int(fields[7]),
                path_end=int(fields[8]),
                path_length=int(fields[6]),
                strand=fields[4],
                mate=int(tags.get("mt", 0)),
                tags=tags,
            )


def write_jsonl(records: Iterator[dict], path: str | Path) -> None:
    with open_text(path, "wt") as handle:
        for record in records:
            handle.write(json.dumps(record, ensure_ascii=False, sort_keys=True) + "\n")


def iter_bam_reads(path: str | Path) -> Iterator[tuple[str, str, int]]:
    """Yield query name, sequence and mate number from BAM/CRAM.

    pysam is an optional dependency. Secondary/supplementary records and records
    without query sequence are skipped. Each physical read is emitted once when
    duplicate alignment records are present.
    """
    try:
        import pysam  # type: ignore
    except ImportError as exc:
        raise RuntimeError("BAM/CRAM input requires: pip install 'str-pg[bio]' or pysam>=0.22") from exc
    seen: set[tuple[str, int]] = set()
    mode = "rc" if str(path).endswith(".cram") else "rb"
    with pysam.AlignmentFile(str(path), mode) as bam:
        for record in bam.fetch(until_eof=True):
            if record.is_secondary or record.is_supplementary or record.query_sequence is None:
                continue
            mate = 1 if record.is_read1 else 2 if record.is_read2 else 0
            key = (record.query_name, mate)
            if key in seen:
                continue
            seen.add(key)
            yield record.query_name, record.query_sequence.upper(), mate
