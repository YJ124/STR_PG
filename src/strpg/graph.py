from __future__ import annotations

import json
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from .io import read_catalog, read_fasta
from .models import AlleleRecord, LocusRecord
from .registry import AlleleRegistry
from .utils import open_text


@dataclass(order=True)
class Event:
    start: int
    end: int
    kind: str
    data: dict


def _parse_candidate_lengths(raw: str, reference_count: int) -> list[int]:
    values = {reference_count}
    for item in (raw or "").replace(";", ",").split(","):
        item = item.strip()
        if item:
            values.add(int(item))
    return sorted(x for x in values if x >= 0)


def _load_simple_vcf(path: str | Path | None) -> dict[str, list[Event]]:
    events: dict[str, list[Event]] = {}
    if path is None:
        return events
    with open_text(path, "rt") as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            chrom, pos_s, vid, ref, alt_s = fields[:5]
            pos = int(pos_s)
            alts = [a for a in alt_s.split(",") if a and not a.startswith("<") and "[" not in a and "]" not in a]
            if not alts:
                continue
            events.setdefault(chrom, []).append(Event(pos, pos + len(ref) - 1, "VAR", {"id": vid or ".", "ref": ref.upper(), "alts": [a.upper() for a in alts]}))
    return events


def build_graph_and_registry(
    reference_path: str | Path,
    catalog_path: str | Path,
    graph_out: str | Path,
    registry_out: str | Path,
    build: str | None = None,
    flank: int = 100,
    vcf_path: str | Path | None = None,
    default_population: str = "GLOBAL",
    initial_count_total: float = 1000.0,
) -> None:
    reference = read_fasta(reference_path)
    rows = read_catalog(catalog_path)
    events_by_chrom = _load_simple_vcf(vcf_path)
    registry = AlleleRegistry()

    for row in rows:
        if build and row["build"] != build:
            continue
        chrom = row["chrom"]
        if chrom not in reference:
            raise ValueError(f"Catalog chromosome {chrom} not in FASTA")
        start, end = int(row["start"]), int(row["end"])
        motif = row["motif"].upper()
        pointer_id = row["pointer_id"] or f"PTR_{chrom}_{start}_{end}"
        if start < 1 or end > len(reference[chrom]) or start > end:
            raise ValueError(f"Invalid coordinates for {pointer_id}")
        ref_repeat = reference[chrom][start - 1 : end]
        ref_count = max(1, round(len(ref_repeat) / len(motif)))
        candidate_lengths = _parse_candidate_lengths(row.get("candidate_Ls", ""), ref_count)
        left = reference[chrom][max(0, start - 1 - flank) : start - 1]
        right = reference[chrom][end : min(len(reference[chrom]), end + flank)]
        freq_payload: dict[str, dict[str, float]] = {}
        if row.get("frequencies_json"):
            freq_payload = json.loads(row["frequencies_json"])
        alleles: list[AlleleRecord] = []
        for repeat_count in candidate_lengths:
            pop_freqs = {
                pop: float(per_allele.get(str(repeat_count), per_allele.get(repeat_count, 0.0)))
                for pop, per_allele in freq_payload.items()
            }
            if not pop_freqs:
                pop_freqs = {default_population: 1.0 / len(candidate_lengths)}
            full_sequence = left + motif * repeat_count + right
            alleles.append(
                AlleleRecord(
                    allele_id=f"{pointer_id}:L{repeat_count}",
                    repeat_count=repeat_count,
                    sequence=full_sequence,
                    repeat_sequence=motif * repeat_count,
                    motif=motif,
                    repeat_start=len(left),
                    repeat_end=len(left) + len(motif) * repeat_count,
                    frequencies=pop_freqs,
                    counts={pop: freq * initial_count_total for pop, freq in pop_freqs.items()},
                    status="registered",
                    source="catalog",
                )
            )
        locus = LocusRecord(
            pointer_id=pointer_id,
            build=row["build"],
            chrom=chrom,
            start=start,
            end=end,
            motif=motif,
            left_flank=left,
            right_flank=right,
            anchor_length=flank,
            alleles=alleles,
            metadata={"reference_repeat_count": ref_count, "reference_repeat_sequence": ref_repeat},
        )
        if pointer_id in registry.loci:
            raise ValueError(f"Duplicate pointer_id: {pointer_id}")
        registry.loci[pointer_id] = locus
        events_by_chrom.setdefault(chrom, []).append(Event(start, end, "STR", {"pointer_id": pointer_id}))

    errors = registry.validate()
    if errors:
        raise ValueError("Registry validation failed:\n" + "\n".join(errors))

    with open_text(graph_out, "wt") as out:
        out.write("H\tVN:Z:1.0\tPG:Z:STR-PG-1.0\n")
        node_id = 0
        for chrom, seq in reference.items():
            events = sorted(events_by_chrom.get(chrom, []))
            filtered: list[Event] = []
            last_end = 0
            for event in events:
                if event.start <= last_end:
                    print(f"[STR-PG] warning: skipping overlapping {event.kind} event {chrom}:{event.start}-{event.end}", file=sys.stderr)
                    continue
                filtered.append(event)
                last_end = event.end

            path_nodes: list[str] = []
            previous_exits: list[str] = []
            cursor = 1

            def emit_segment(sequence: str, tags: list[str] | None = None) -> str:
                nonlocal node_id
                node_id += 1
                nid = f"n{node_id}"
                seq_field = sequence if sequence else "N"
                out.write("\t".join(["S", nid, seq_field] + (tags or [])) + "\n")
                return nid

            def connect(srcs: Iterable[str], dsts: Iterable[str]) -> None:
                for src in srcs:
                    for dst in dsts:
                        out.write(f"L\t{src}\t+\t{dst}\t+\t0M\n")

            for event in filtered:
                if event.start > cursor:
                    chunk = emit_segment(seq[cursor - 1 : event.start - 1], [f"SN:Z:{chrom}", f"SO:i:{cursor-1}"])
                    if previous_exits:
                        connect(previous_exits, [chunk])
                    path_nodes.append(chunk)
                    previous_exits = [chunk]
                if event.kind == "STR":
                    pid = event.data["pointer_id"]
                    ptr = emit_segment(
                        "N",
                        ["TP:Z:STR_POINTER", f"PI:Z:{pid}", f"SN:Z:{chrom}", f"SO:i:{event.start-1}", f"EO:i:{event.end}"],
                    )
                    if previous_exits:
                        connect(previous_exits, [ptr])
                    path_nodes.append(ptr)
                    previous_exits = [ptr]
                else:
                    ref_node = emit_segment(event.data["ref"], ["TP:Z:VAR_REF", f"SN:Z:{chrom}", f"SO:i:{event.start-1}"])
                    alt_nodes = [emit_segment(alt, ["TP:Z:VAR_ALT", f"VI:Z:{event.data['id']}"]) for alt in event.data["alts"]]
                    if previous_exits:
                        connect(previous_exits, [ref_node] + alt_nodes)
                    path_nodes.append(ref_node)
                    previous_exits = [ref_node] + alt_nodes
                cursor = event.end + 1

            if cursor <= len(seq):
                tail = emit_segment(seq[cursor - 1 :], [f"SN:Z:{chrom}", f"SO:i:{cursor-1}"])
                if previous_exits:
                    connect(previous_exits, [tail])
                path_nodes.append(tail)
                previous_exits = [tail]
            if path_nodes:
                out.write(f"P\t{chrom}\t{','.join(n + '+' for n in path_nodes)}\t*\n")
        for pid, locus in sorted(registry.loci.items()):
            out.write(
                "X\tPTR\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                    pid, locus.chrom, locus.start, locus.end, locus.motif, locus.anchor_length
                )
            )
    registry.save(registry_out)
