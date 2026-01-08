#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PGG (Pan-Genome Graph) — Build step (4.1):

输入：
  - 参考基因组 FASTA（--ref）
  - 变异 VCF（--vcf，可选；支持 SNP/小 indel；与 STR 不重叠为宜）
  - STR 位点目录 TSV（--str-catalog，可选；初期建议用 TSV 提供 STR 信息）

输出：
  - GFA-like 图文件（--out graph.gfa），包含：
      S（节点）、L（边）、P（路径）、X（扩展元数据：X STR / X PTR）

要点：
  - 每条染色体生成一个“参考路径 P <chrom> ...”，沿参考序列穿过所有变异位点；
  - 每个变异在图上形成“泡泡”（bubble）：左锚节点 →（等位路径）→ 右锚/间隔节点；
  - STR：为每个候选重复次数 L 生成一条等位路径，路径名含 L，写入 `X PTR <pid> <pointer_id>`；
  - 为每个 STR 位点写 `X STR <locus_id> <motif> <left_flank> <right_flank> <candidate_Ls>`；

依赖：
  - Python 3.9+
  - pysam（解析 VCF；若不用 VCF，可不安装）

示例：
  python pgg_build.py \
    --ref ref.fa \
    --vcf variants.vcf.gz \
    --str-catalog str_catalog.tsv \
    --build GRCh38 \
    --flank 100 \
    --out graph.gfa

STR TSV 期望列（制表符分隔）：
  build  chrom  start  end  motif  candidate_Ls  pointer_id
  - build: 参考版本（如 GRCh38）
  - start/end：1-based 闭区间
  - motif：如 "CA"
  - candidate_Ls：逗号分隔的整数，如 "8,9,10,11,13"（可留空，自动含参考 L）
  - pointer_id：可留空，将自动生成

注意：
  - 初版不处理变异间复杂重叠；若发现重叠会跳过后者并给出告警。
  - GFA L 行的 overlap 字段统一写 "0M"。
"""
from __future__ import annotations
import argparse
import sys
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional, Iterable, Set

# 可选依赖：仅在 --vcf 提供时需要
try:
    import pysam  # type: ignore
except Exception:
    pysam = None

# =========================
# 基础数据结构
# =========================
@dataclass
class Node:
    id: int
    seq: str  # 序列（A/C/G/T/N）；不建议零长

@dataclass
class Edge:
    u: int
    v: int

@dataclass
class Event:
    """一个变异事件（SNP/INDEL 或 STR）"""
    chrom: str
    start: int  # 1-based inclusive
    end: int    # 1-based inclusive
    kind: str   # 'snpindel' | 'str'
    
    # SNP/indel 专用
    ref: Optional[str] = None
    alts: Optional[List[str]] = None  # ALT 序列列表（多等位时多条）

    # STR 专用
    motif: Optional[str] = None
    candidate_Ls: Optional[List[int]] = None
    pointer_id: Optional[str] = None  # 指向频率库的位点主键

# =========================
# FASTA I/O（无依赖）
# =========================

def read_fasta(path: str) -> Dict[str, str]:
    """读取（小规模）FASTA 为内存字典：{chrom: seq}
    注意：大参考建议用 faidx/随机访问；此处为原型实现。
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
# STR 目录 TSV 读取
# =========================

def parse_str_catalog(tsv_path: str) -> Dict[Tuple[str,int,int], dict]:
    """读取 STR 目录：返回 {(chrom,start,end): info_dict}
    需要列：build, chrom, start, end, motif, candidate_Ls(可空), pointer_id(可空)
    """
    import csv
    out: Dict[Tuple[str,int,int], dict] = {}
    with open(tsv_path, 'r') as f:
        rd = csv.DictReader(f, delimiter='\t')
        required = {'build','chrom','start','end','motif'}
        miss = required - set(rd.fieldnames or [])
        if miss:
            raise ValueError(f"STR catalog缺少列: {sorted(miss)}")
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
                    raise ValueError(f"candidate_Ls 解析失败: {cand_raw}")
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
# VCF 解析（SNP/indel）
# =========================

def parse_vcf_events(vcf_path: str) -> Dict[str, List[Event]]:
    """读取 VCF，抽取 SNP/小 indel 事件（不含 STR 特殊逻辑）。
    仅当安装了 pysam 才可用。
    返回：按染色体分组的 Event 列表，已按 start 排序。
    """
    if vcf_path is None:
        return {}
    if pysam is None:
        raise RuntimeError("需要安装 pysam 才能解析 VCF。pip install pysam")
    out: Dict[str, List[Event]] = {}
    vf = pysam.VariantFile(vcf_path)
    for rec in vf.fetch():  # 需有 .tbi 或 .csi 索引；未索引则遍历
        chrom = str(rec.chrom)
        pos = int(rec.pos)  # 1-based
        ref = str(rec.ref).upper()
        alts = [str(a).upper() for a in (rec.alts or []) if a and a != '<*>' and a != '*']
        if not alts:
            continue
        # 简化：跳过过长 ALT（例如 SV），限制到小 indel
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
    # 排序
    for c in out:
        out[c].sort(key=lambda e: e.start)
    return out

# =========================
# 构图器
# =========================
class GraphBuilder:
    def __init__(self, build: str, flank: int = 100):
        self.build = build
        self.flank = flank
        self.nodes: List[Node] = []
        self.edges: Set[Tuple[int,int]] = set()
        self.paths: List[Tuple[str, List[int]]] = []  # (path_id, [node_ids])
        self.x_lines: List[str] = []  # 原样写入的 X 行
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
        """在单条染色体上构图：
        - 参考路径：穿越各事件的 REF 分支
        - 每个事件：构造泡泡；STR 生成多条 L 分支与 X 行
        """
        # 按位置排序 + 去重/跳过重叠
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

        # 预处理：若有事件则先建“头部”左锚节点
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
            # 无事件：整条序列一个节点即可
            nid = self.add_node(ref_seq)
            self.add_path(chrom, [nid])
            return

        # 逐事件生成泡泡
        for i, ev in enumerate(evs):
            # 计算右侧间隔（到下一个事件前）
            next_start = evs[i+1].start if i+1 < len(evs) else (chrom_len + 1)
            gap_left = ev.end + 1
            gap_right = next_start - 1
            gap_seq = ref_seq[gap_left-1:gap_right] if gap_right >= gap_left else ''

            # —— REF 分支节点 ——
            ref_seg_seq = ref_seq[ev.start-1:ev.end]
            ref_node = self.add_node(ref_seg_seq)
            # 左锚 → REF
            self.add_edge(left_anchor, ref_node)

            # —— ALT/等位分支节点 ——
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
                    # 为该 ALT 生成一个局部路径（左锚 → ALT → GAP）稍后补 GAP
            elif ev.kind == 'str':
                assert ev.motif is not None
                motif = ev.motif.upper()
                # 参考 L
                k = len(motif)
                L_ref = len(ref_seg_seq) // k if k>0 else 0
                # 候选 L 集合：目录给的 + 参考 L
                cand = set(ev.candidate_Ls or [])
                if L_ref > 0:
                    cand.add(L_ref)
                candidate_Ls = sorted(cand)

                # X STR 行（位点注册）
                left_flank = ref_seq[max(0, ev.start-1-self.flank):ev.start-1]
                right_flank = ref_seq[ev.end: min(chrom_len, ev.end+self.flank)]
                locus_id = self._locus_id(self.build, chrom, ev.start, ev.end, motif)
                x_str = (
                    f"X\tSTR\t{locus_id}\t{motif}\t{left_flank}\t{right_flank}\t"
                    + ",".join(str(x) for x in candidate_Ls)
                )
                self.x_lines.append(x_str)

                # pointer_id（频率库主键）
                pointer_id = ev.pointer_id or locus_id

                for L in candidate_Ls:
                    alt_seq = motif * L
                    if alt_seq == ref_seg_seq:
                        continue  # 与 REF 相同无需另一分支
                    alt_node = self.add_node(alt_seq)
                    self.add_edge(left_anchor, alt_node)
                    pid = f"allele_{chrom}_{ev.start}_{ev.end}_{motif}_L{L}"
                    alt_nodes.append((pid, alt_node))
                    # X PTR：将此等位路径 pid 绑定 pointer_id
                    self.x_lines.append(f"X\tPTR\t{pid}\t{pointer_id}")
            else:
                raise ValueError(f"unknown event kind: {ev.kind}")

            # —— 右侧 GAP 节点 ——
            gap_node = None
            if gap_seq:
                gap_node = self.add_node(gap_seq)
            # 连接 REF / ALT → GAP
            self.add_edge(ref_node, gap_node)
            for pid, aid in alt_nodes:
                self.add_edge(aid, gap_node)
                # 把 ALT 的局部路径记录出来：左锚 → ALT → GAP
                local_path_nodes = []
                if left_anchor is not None:
                    local_path_nodes.append(left_anchor)
                local_path_nodes.append(aid)
                if gap_node is not None:
                    local_path_nodes.append(gap_node)
                self.add_path(pid, local_path_nodes)

            # 参考路径：追加 REF → GAP
            ref_path_nodes.append(ref_node)
            if gap_node is not None:
                ref_path_nodes.append(gap_node)
            # 下一个事件的左锚 = 当前 GAP（没有 GAP 则为 None，表示紧邻）
            left_anchor = gap_node

        # 写整条染色体的参考路径
        self.add_path(chrom, ref_path_nodes)

    # ----
    def write_gfa(self, out_path: str):
        with open(out_path, 'w') as w:
            # Header（可选）：写个版本标记
            w.write("H\tVN:Z:1.0\tpg:Z:pgg-build\n")
            # S 行
            for n in self.nodes:
                seq = n.seq if n.seq else "*"
                w.write(f"S\t{n.id}\t{seq}\n")
            # L 行（去重后写）
            for u, v in sorted(self.edges):
                if v is None or u is None:
                    continue
                w.write(f"L\t{u}\t+\t{v}\t+\t0M\n")
            # P 行
            for pid, nodes in self.paths:
                if not nodes:
                    continue
                nds = ",".join(f"{nid}+" for nid in nodes)
                w.write(f"P\t{pid}\t{nds}\t*\n")
            # X 行（扩展）
            for line in self.x_lines:
                w.write(line + "\n")

# =========================
# 事件汇总（合并 VCF 与 STR 目录）
# =========================

def combine_events_for_chrom(chrom: str,
                             vcf_events: Dict[str,List[Event]],
                             str_catalog: Dict[Tuple[str,int,int],dict]) -> List[Event]:
    evs: List[Event] = []
    # VCF（SNP/indel）
    for e in vcf_events.get(chrom, []):
        evs.append(e)
    # STR 目录
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
# 主流程
# =========================

def main():
    ap = argparse.ArgumentParser(description='PGG build: 从 FASTA + VCF/STR 目录构建 GFA 图')
    ap.add_argument('--ref', required=True, help='参考 FASTA')
    ap.add_argument('--vcf', default=None, help='变异 VCF（可选）')
    ap.add_argument('--str-catalog', default=None, help='STR 目录 TSV（可选，推荐）')
    ap.add_argument('--build', required=True, help='参考版本名（用于 pointer_id / locus_id，如 GRCh38）')
    ap.add_argument('--flank', type=int, default=100, help='X STR 左右侧翼长度（默认100）')
    ap.add_argument('--out', required=True, help='输出 GFA 文件')
    args = ap.parse_args()

    # 读 FASTA
    sys.stderr.write('[INFO] reading FASTA...\n')
    ref = read_fasta(args.ref)
    if not ref:
        raise SystemExit('FASTA 为空。')

    # 读 VCF 事件
    if args.vcf:
        sys.stderr.write('[INFO] parsing VCF events...\n')
        vcf_events = parse_vcf_events(args.vcf)
    else:
        vcf_events = {}

    # 读 STR 目录
    if args.str_catalog:
        sys.stderr.write('[INFO] parsing STR catalog...\n')
        str_catalog = parse_str_catalog(args.str_catalog)
    else:
        str_catalog = {}

    # 构图
    gb = GraphBuilder(build=args.build, flank=args.flank)

    # 每条染色体独立处理
    for chrom, seq in ref.items():
        sys.stderr.write(f'[INFO] building chrom {chrom}...\n')
        evs = combine_events_for_chrom(chrom, vcf_events, str_catalog)
        gb.build_chrom(chrom, seq, evs, str_catalog)

    # 输出 GFA
    sys.stderr.write(f'[INFO] writing GFA to {args.out}\n')
    gb.write_gfa(args.out)
    sys.stderr.write('[INFO] done.\n')

if __name__ == '__main__':
    main()
