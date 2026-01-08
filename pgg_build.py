#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PGG (Pan-Genome Graph) — Build step (4.1):

Author：YJ

Inputs:
  - Reference Genome FASTA (--ref)
  - Variant VCF (--vcf, optional; supports SNP/small indels; preferably non-overlapping with STRs)
  - STR Locus Catalog TSV (--str-catalog, optional; recommended for providing STR info initially)

Outputs:
  - GFA-like graph file (--out graph.gfa), containing:
      S (Segments/Nodes), L (Links/Edges), P (Paths), X (Extended metadata: X STR / X PTR)

Key Points:
  - A "Reference Path P <chrom> ..." is generated for each chromosome, traversing the reference sequence through all variant sites.
  - Each variant forms a "bubble" on the graph: Left Anchor Node → (Allele Paths) → Right Anchor/Gap Node.
  - STR: Generates an allele path for each candidate repeat count L. Path name includes L. Writes `X PTR <pid> <pointer_id>`.
  - Writes `X STR <locus_id> <motif> <left_flank> <right_flank> <candidate_Ls>` for each STR locus.

Dependencies:
  - Python 3.9+
  - pysam (for parsing VCF; not required if VCF is not used)

Example:
  python pgg_build.py \
    --ref ref.fa \
    --vcf variants.vcf.gz \
    --str-catalog str_catalog.tsv \
    --build GRCh38 \
    --flank 100 \
    --out graph.gfa

STR TSV Expected Columns (Tab-separated):
  build  chrom  start  end  motif  candidate_Ls  pointer_id
  - build: Reference version (e.g., GRCh38)
  - start/end: 1-based inclusive
  - motif: e.g., "CA"
  - candidate_Ls: Comma-separated integers, e.g., "8,9,10,11,13" (can be empty, automatically includes reference L)
  - pointer_id: Can be empty, will be automatically generated

Notes:
  - The initial version does not handle complex overlaps between variants; if overlap is found, the latter is skipped with a warning.
  - The overlap field in GFA L lines is uniformly set to "0M".
"""
from __future__ import annotations
import argparse
import sys
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional, Iterable, Set

# Optional dependency: required only when --vcf is provided
try:
    import pysam  # type: ignore
except Exception:
    pysam = None

# =========================
# Basic Data Structures
# =========================
@dataclass
class Node:
    id: int
    seq: str  # Sequence (A/C/G/T/N); zero-length not recommended

@dataclass
class Edge:
    u: int
    v: int

@dataclass
class Event:
    """A variant event (SNP/INDEL or STR)"""
    chrom: str
    start: int  # 1-based inclusive
    end: int    # 1-based inclusive
    kind: str   # 'snpindel' | 'str'
    
    # Specific to SNP/indel
    ref: Optional[str] = None
    alts: Optional[List[str]] = None  # List of ALT sequences (multiple for multi-allelic)

    # Specific to STR
    motif: Optional[str] = None
    candidate_Ls: Optional[List[int]] = None
    pointer_id: Optional[str] = None  # Locus primary key pointing to frequency db

# =========================
# FASTA I/O (No dependencies)
# =========================

def read_fasta(path: str) -> Dict[str, str]:
    """Read (small-scale) FASTA into memory dict: {chrom: seq}
    Note: For large references, use faidx/random access; this is a prototype implementation.
    """
    seqs: Dict[str, List[str]] = {}
    name: Optional[str] = None
    with open(path, 'r') as f:
        for line in f:
            if not line:
                continue
            if line.startswith('>'):
                name = line[1:].strip().split()[0]
                if name in seqs:
                    raise ValueError(f"FASTA has duplicate contig: {name}")
                seqs[name] = []
            else:
                if name is None:
                    raise ValueError("FASTA format error: sequence before header")
                seqs[name].append(line.strip().upper())
    return {k: ''.join(v) for k, v in seqs.items()}

# =========================
# STR Catalog TSV Reading
# =========================

def parse_str_catalog(tsv_path: str) -> Dict[Tuple[str,int,int], dict]:
    """Read STR catalog: returns {(chrom,start,end): info_dict}
    Required columns: build, chrom, start, end, motif, candidate_Ls(nullable), pointer_id(nullable)
    """
    import csv
    out: Dict[Tuple[str,int,int], dict] = {}
    with open(tsv_path, 'r') as f:
        rd = csv.DictReader(f, delimiter='\t')
        required = {'build','chrom','start','end','motif'}
        miss = required - set(rd.fieldnames or [])
        if miss:
            raise ValueError(f"STR catalog missing columns: {sorted(miss)}")
        for row in rd:
            chrom = row['chrom']
            start = int(row['start'])
            end = int(row['end'])
            motif = row['motif'].upper()
            cand_raw = (row.get('candidate_Ls') or '').strip()
            candidate_Ls = None
            if cand_raw:
                try:
                    candidate_Ls = [int(x) for x in cand_raw.replace(' ','').split(',') if x]
                except Exception:
                    raise ValueError(f"Failed to parse candidate_Ls: {cand_raw}")
            pointer_id = (row.get('pointer_id') or '').strip() or None
            build = (row.get('build') or '').strip()
            out[(chrom,start,end)] = {
                'build': build,
                'chrom': chrom,
                'start': start,
                'end': end,
                'motif': motif,
                'candidate_Ls': candidate_Ls,
                'pointer_id': pointer_id,
            }
    return out

# =========================
# VCF Parsing (SNP/indel)
# =========================

def parse_vcf_events(vcf_path: str) -> Dict[str, List[Event]]:
    """Read VCF, extract SNP/small indel events (excluding STR special logic).
    Only available if pysam is installed.
    Returns: List of Events grouped by chromosome, sorted by start.
    """
    if vcf_path is None:
        return {}
    if pysam is None:
        raise RuntimeError("pysam is required to parse VCF. pip install pysam")
    out: Dict[str, List[Event]] = {}
    vf = pysam.VariantFile(vcf_path)
    for rec in vf.fetch():  # Requires .tbi or .csi index; iterates if not indexed
        chrom = str(rec.chrom)
        pos = int(rec.pos)  # 1-based
        ref = str(rec.ref).upper()
        alts = [str(a).upper() for a in (rec.alts or []) if a and a != '<*>' and a != '*']
        if not alts:
            continue
        # Simplification: Skip overly long ALTs (e.g., SVs), restrict to small indels
        max_len = max(len(ref), *(len(a) for a in alts))
        if max_len > 50:
            continue
        ev = Event(
            chrom=chrom,
            start=pos,
            end=pos + len(ref) - 1,
            kind='snpindel',
            ref=ref,
            alts=alts,
        )
        out.setdefault(chrom, []).append(ev)
    # Sort
    for c in out:
        out[c].sort(key=lambda e: e.start)
    return out

# =========================
# Graph Builder
# =========================
class GraphBuilder:
    def __init__(self, build: str, flank: int = 100):
        self.build = build
        self.flank = flank
        self.nodes: List[Node] = []
        self.edges: Set[Tuple[int,int]] = set()
        self.paths: List[Tuple[str, List[int]]] = []  # (path_id, [node_ids])
        self.x_lines: List[str] = []  # X lines written as-is
        self._nid = 0

    def add_node(self, seq: str) -> int:
        assert len(seq) >= 0
        self._nid += 1
        nid = self._nid
        self.nodes.append(Node(nid, seq))
        return nid

    def add_edge(self, u: Optional[int], v: Optional[int]):
        if u is None or v is None:
            return
        if u == v:
            return
        self.edges.add((u, v))

    def add_path(self, pid: str, node_ids: List[int]):
        self.paths.append((pid, node_ids))

    # ----
    @staticmethod
    def _locus_id(build: str, chrom: str, start: int, end: int, motif: str) -> str:
        return f"{build}:{chrom}:{start}-{end}:{motif}"

    def build_chrom(self,
                    chrom: str,
                    ref_seq: str,
                    events: List[Event],
                    str_catalog: Dict[Tuple[str,int,int], dict]):
        """Build graph on a single chromosome:
        - Reference Path: REF branch traversing each event
        - Each Event: Construct bubble; STR generates multiple L branches and X lines
        """
        # Sort by position + deduplicate/skip overlaps
        evs = sorted(events, key=lambda e: e.start) if events else []
        filtered: List[Event] = []
        last_end = 0
        for ev in evs:
            if ev.start <= last_end:
                sys.stderr.write(f"[WARN] overlap detected at {chrom}:{ev.start}-{ev.end}, skip.\n")
                continue
            filtered.append(ev)
            last_end = ev.end
        evs = filtered

        chrom_len = len(ref_seq)
        ref_path_nodes: List[int] = []

        # Preprocessing: If events exist, build the "head" left anchor node first
        left_anchor: Optional[int] = None
        head_start = 1
        if evs:
            first = evs[0]
            if first.start > 1:
                head_seq = ref_seq[0:first.start-1]
                if head_seq:
                    left_anchor = self.add_node(head_seq)
                    ref_path_nodes.append(left_anchor)
        else:
            # No events: The entire sequence is just one node
            nid = self.add_node(ref_seq)
            self.add_path(chrom, [nid])
            return

        # Generate bubbles event by event
        for i, ev in enumerate(evs):
            # Calculate right gap (until next event)
            next_start = evs[i+1].start if i+1 < len(evs) else (chrom_len + 1)
            gap_left = ev.end + 1
            gap_right = next_start - 1
            gap_seq = ref_seq[gap_left-1:gap_right] if gap_right >= gap_left else ''

            # —— REF branch node ——
            ref_seg_seq = ref_seq[ev.start-1:ev.end]
            ref_node = self.add_node(ref_seg_seq)
            # Left Anchor → REF
            self.add_edge(left_anchor, ref_node)

            # —— ALT/Allele branch nodes ——
            alt_nodes: List[Tuple[str,int]] = []  # (path_id, node_id)

            if ev.kind == 'snpindel':
                assert ev.alts is not None
                for ai, alt in enumerate(ev.alts):
                    if alt == ref_seg_seq:
                        continue
                    alt_node = self.add_node(alt)
                    self.add_edge(left_anchor, alt_node)
                    pid = f"alt_{chrom}_{ev.start}_{ev.end}_a{ai+1}"
                    alt_nodes.append((pid, alt_node))
                    # Generate a local path for this ALT (Left Anchor → ALT → GAP), GAP added later
            elif ev.kind == 'str':
                assert ev.motif is not None
                motif = ev.motif.upper()
                # Reference L
                k = len(motif)
                L_ref = len(ref_seg_seq) // k if k>0 else 0
                # Candidate L set: Catalog provided + Reference L
                cand = set(ev.candidate_Ls or [])
                if L_ref > 0:
                    cand.add(L_ref)
                candidate_Ls = sorted(cand)

                # X STR line (Locus registration)
                left_flank = ref_seq[max(0, ev.start-1-self.flank):ev.start-1]
                right_flank = ref_seq[ev.end: min(chrom_len, ev.end+self.flank)]
                locus_id = self._locus_id(self.build, chrom, ev.start, ev.end, motif)
                x_str = (
                    f"X\tSTR\t{locus_id}\t{motif}\t{left_flank}\t{right_flank}\t"
                    + ",".join(str(x) for x in candidate_Ls)
                )
                self.x_lines.append(x_str)

                # pointer_id (Primary key for frequency db)
                pointer_id = ev.pointer_id or locus_id

                for L in candidate_Ls:
                    alt_seq = motif * L
                    if alt_seq == ref_seg_seq:
                        continue  # No need for another branch if identical to REF
                    alt_node = self.add_node(alt_seq)
                    self.add_edge(left_anchor, alt_node)
                    pid = f"allele_{chrom}_{ev.start}_{ev.end}_{motif}_L{L}"
                    alt_nodes.append((pid, alt_node))
                    # X PTR: Bind this allele path pid to pointer_id
                    self.x_lines.append(f"X\tPTR\t{pid}\t{pointer_id}")
            else:
                raise ValueError(f"unknown event kind: {ev.kind}")

            # —— Right side GAP node ——
            gap_node = None
            if gap_seq:
                gap_node = self.add_node(gap_seq)
            # Connect REF / ALT → GAP
            self.add_edge(ref_node, gap_node)
            for pid, aid in alt_nodes:
                self.add_edge(aid, gap_node)
                # Record the local path for ALT: Left Anchor → ALT → GAP
                local_path_nodes = []
                if left_anchor is not None:
                    local_path_nodes.append(left_anchor)
                local_path_nodes.append(aid)
                if gap_node is not None:
                    local_path_nodes.append(gap_node)
                self.add_path(pid, local_path_nodes)

            # Reference Path: Append REF → GAP
            ref_path_nodes.append(ref_node)
            if gap_node is not None:
                ref_path_nodes.append(gap_node)
            # Next event's left anchor = Current GAP (None if no GAP, meaning adjacent)
            left_anchor = gap_node

        # Write the reference path for the whole chromosome
        self.add_path(chrom, ref_path_nodes)

    # ----
    def write_gfa(self, out_path: str):
        with open(out_path, 'w') as w:
            # Header (Optional): Write a version tag
            w.write("H\tVN:Z:1.0\tpg:Z:pgg-build\n")
            # S lines
            for n in self.nodes:
                seq = n.seq if n.seq else "*"
                w.write(f"S\t{n.id}\t{seq}\n")
            # L lines (Write after deduplication)
            for u, v in sorted(self.edges):
                if v is None or u is None:
                    continue
                w.write(f"L\t{u}\t+\t{v}\t+\t0M\n")
            # P lines
            for pid, nodes in self.paths:
                if not nodes:
                    continue
                nds = ",".join(f"{nid}+" for nid in nodes)
                w.write(f"P\t{pid}\t{nds}\t*\n")
            # X lines (Extended)
            for line in self.x_lines:
                w.write(line + "\n")

# =========================
# Event Summarization (Merge VCF and STR Catalog)
# =========================

def combine_events_for_chrom(chrom: str,
                             vcf_events: Dict[str,List[Event]],
                             str_catalog: Dict[Tuple[str,int,int],dict]) -> List[Event]:
    evs: List[Event] = []
    # VCF (SNP/indel)
    for e in vcf_events.get(chrom, []):
        evs.append(e)
    # STR Catalog
    for (c,s,e), info in str_catalog.items():
        if c != chrom:
            continue
        evs.append(Event(
            chrom=chrom,
            start=s,
            end=e,
            kind='str',
            motif=info['motif'],
            candidate_Ls=info.get('candidate_Ls'),
            pointer_id=info.get('pointer_id') or GraphBuilder._locus_id(info.get('build') or '', chrom, s, e, info['motif'])
        ))
    evs.sort(key=lambda x: x.start)
    return evs

# =========================
# Main Workflow
# =========================

def main():
    ap = argparse.ArgumentParser(description='PGG build: Build GFA graph from FASTA + VCF/STR Catalog')
    ap.add_argument('--ref', required=True, help='Reference FASTA')
    ap.add_argument('--vcf', default=None, help='Variant VCF (optional)')
    ap.add_argument('--str-catalog', default=None, help='STR Catalog TSV (optional, recommended)')
    ap.add_argument('--build', required=True, help='Reference build name (for pointer_id / locus_id, e.g., GRCh38)')
    ap.add_argument('--flank', type=int, default=100, help='X STR left/right flank length (default 100)')
    ap.add_argument('--out', required=True, help='Output GFA file')
    args = ap.parse_args()

    # Read FASTA
    sys.stderr.write('[INFO] reading FASTA...\n')
    ref = read_fasta(args.ref)
    if not ref:
        raise SystemExit('FASTA is empty.')

    # Read VCF events
    if args.vcf:
        sys.stderr.write('[INFO] parsing VCF events...\n')
        vcf_events = parse_vcf_events(args.vcf)
    else:
        vcf_events = {}

    # Read STR catalog
    if args.str_catalog:
        sys.stderr.write('[INFO] parsing STR catalog...\n')
        str_catalog = parse_str_catalog(args.str_catalog)
    else:
        str_catalog = {}

    # Build Graph
    gb = GraphBuilder(build=args.build, flank=args.flank)

    # Process each chromosome independently
    for chrom, seq in ref.items():
        sys.stderr.write(f'[INFO] building chrom {chrom}...\n')
        evs = combine_events_for_chrom(chrom, vcf_events, str_catalog)
        gb.build_chrom(chrom, seq, evs, str_catalog)

    # Output GFA
    sys.stderr.write(f'[INFO] writing GFA to {args.out}\n')
    gb.write_gfa(args.out)
    sys.stderr.write('[INFO] done.\n')

if __name__ == '__main__':
    main()
