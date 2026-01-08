#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PGG（Pan-Genome Graph）— 4.2 index（无 SQL 版）

将 build 产出的 GFA 建立以下索引：
  1) 路径扁平化序列池（每条 P 路径一条 .fa）
  2) 路径坐标索引：{pid -> [(start,node_id,len), ...]}（pickle）
  3) syncmer/minimizer 种子索引：k-mer -> postings(path_id, pos)
     · 采用 Python 标准库 dbm 的分片 KV 存储（无 SQL，无第三方依赖）
     · key=8字节hash（blake2b），value=pickle(postings 列表)
  4) 位点倒排：locus_id -> {allele_paths, candidate_Ls, anchor区段}

用法：
  python pgg_index.py \
    --graph graph.gfa \
    --out index_dir \
    --method syncmer --k 15 --s 5 --t 2 \
    --shards 64 --anchor 100

输出目录结构：
  index_dir/
    ├─ paths/                 # 扁平化序列池
    ├─ path_coord.pkl         # 路径坐标索引
    ├─ locus.json             # 位点倒排索引
    ├─ pid_map.json           # 路径名 <-> 整数ID 映射
    └─ seeds/                 # 分片KV（dbm）
        ├─ seeds_000.dbm      # （实际可能含 .db/.dat/.dir 等，取决于平台）
        ├─ seeds_001.dbm
        └─ meta.json          # 方法与参数、路径长度、shard 数等
"""
from __future__ import annotations
import os, sys, re, argparse, json, pickle, hashlib, collections, struct
from dataclasses import dataclass
from typing import Dict, List, Tuple

# ------------------------
# 通用工具
# ------------------------
DNA_COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')

def revcomp(seq: str) -> str:
    return seq.translate(DNA_COMP)[::-1]

# ------------------------
# 简化的 GFA 解析（兼容 pgg_build.py 输出）
# ------------------------
@dataclass
class Node:
    id: int
    seq: str

class GFA:
    def __init__(self, path: str):
        self.nodes: Dict[int, Node] = {}
        self.paths: Dict[str, List[int]] = {}
        self.x_ptr: Dict[str, str] = {}  # pid -> pointer_id(locus_id)
        self.x_str: Dict[str, dict] = {} # locus_id -> info
        self._parse(path)

    def _parse(self, path: str):
        with open(path, 'r') as f:
            for line in f:
                if not line.strip():
                    continue
                tag = line[0]
                if tag == 'S':
                    parts = line.rstrip().split('\t')
                    nid = int(parts[1]); seq = parts[2]
                    self.nodes[nid] = Node(nid, '' if seq=='*' else seq)
                elif tag == 'P':
                    parts = line.rstrip().split('\t')
                    pid = parts[1]
                    ids = [int(x[:-1]) for x in parts[2].split(',')]
                    self.paths[pid] = ids
                elif line.startswith('X\tPTR'):
                    _, _, pid, pointer = line.rstrip().split('\t')
                    self.x_ptr[pid] = pointer
                elif line.startswith('X\tSTR'):
                    parts = line.rstrip().split('\t')
                    locus_id, motif, lf, rf, cands = parts[2:7]
                    cL = [int(x) for x in cands.split(',')] if cands else []
                    self.x_str[locus_id] = {
                        'motif': motif, 'left_flank': lf, 'right_flank': rf, 'candidate_Ls': cL
                    }

# ------------------------
# 路径扁平化与坐标映射
# ------------------------
@dataclass
class PathSeq:
    pid: str
    seq: str
    segments: List[Tuple[int,int,int]]  # (start, node_id, length)

class PathFlattener:
    @staticmethod
    def flatten(gfa: GFA) -> Dict[str, PathSeq]:
        out: Dict[str, PathSeq] = {}
        for pid, node_ids in gfa.paths.items():
            pos = 0; segs: List[Tuple[int,int,int]] = []; buf = []
            for nid in node_ids:
                s = gfa.nodes[nid].seq
                if s:
                    segs.append((pos, nid, len(s)))
                    buf.append(s); pos += len(s)
            out[pid] = PathSeq(pid, ''.join(buf), segs)
        return out

# ------------------------
# syncmer / minimizer 生成
# ------------------------
class Seeder:
    @staticmethod
    def canonical(kmer: str) -> str:
        rc = revcomp(kmer)
        return kmer if kmer <= rc else rc

    @staticmethod
    def kmer_hash(kmer: str) -> int:
        h = hashlib.blake2b(kmer.encode('ascii'), digest_size=8).digest()
        return int.from_bytes(h, 'big', signed=False)

    @staticmethod
    def syncmer_positions(seq: str, k: int, s: int, t: int):
        n = len(seq)
        if k>n or s>k or t<0 or t>k-s: return []
        pos = []
        for i in range(0, n-k+1):
            kmer = seq[i:i+k]
            minv = None; minp = None
            for j in range(0, k-s+1):
                v = kmer[j:j+s]
                if (minv is None) or (v < minv):
                    minv, minp = v, j
            if minp == t:
                pos.append(i)
        return pos

    @staticmethod
    def minimizer_positions(seq: str, k: int, w: int):
        n = len(seq)
        if k>n or w<=0: return []
        pos = []
        window = []  # (hash, i)
        def kmer(i):
            return Seeder.canonical(seq[i:i+k])
        for i in range(0, n-k+1):
            h = Seeder.kmer_hash(kmer(i))
            window.append((h,i))
            left = i - w + 1
            while window and window[0][1] < left:
                window.pop(0)
            if left >= 0:
                m = min(window, key=lambda x: x[0])
                if not pos or pos[-1] != m[1]:
                    pos.append(m[1])
        return pos

# ------------------------
# DBM 分片 KV（无 SQL）
# ------------------------
class DBMShardWriter:
    def __init__(self, dir_path: str, shards: int):
        import dbm
        self.dir = dir_path
        self.shards = shards
        os.makedirs(self.dir, exist_ok=True)
        self.buffers = [collections.defaultdict(list) for _ in range(shards)]
        self.meta = {
            'version': 1,
            'shards': shards,
            'paths': {},   # pid -> length
            'pid_to_i': {},
            'i_to_pid': {},
            'method': None, 'k': None, 's': None, 't': None, 'w': None,
        }

    def add_path_meta(self, pid: str, length: int, pid_id: int):
        self.meta['paths'][pid] = length
        self.meta['pid_to_i'][pid] = pid_id
        self.meta['i_to_pid'][str(pid_id)] = pid

    def _shard_id(self, key_hash: int) -> int:
        return key_hash % self.shards

    @staticmethod
    def _key_bytes(key_hash: int) -> bytes:
        return struct.pack('>Q', key_hash)  # 8-byte BE

    def add_seed(self, key_hash: int, pid_i: int, pos: int):
        sid = self._shard_id(key_hash)
        self.buffers[sid][key_hash].append((pid_i, pos))

    def finalize(self):
        import dbm
        # 写 meta
        json.dump(self.meta, open(os.path.join(self.dir, 'meta.json'), 'w'))
        # 逐 shard 落盘
        for sid in range(self.shards):
            shard_path = os.path.join(self.dir, f'seeds_{sid:03d}.dbm')
            db = dbm.open(shard_path, 'n')  # new
            try:
                for kh, postings in self.buffers[sid].items():
                    db[self._key_bytes(kh)] = pickle.dumps(postings, protocol=pickle.HIGHEST_PROTOCOL)
            finally:
                db.close()

class DBMShardReader:
    def __init__(self, dir_path: str):
        import dbm
        self.dir = dir_path
        self.meta = json.load(open(os.path.join(self.dir, 'meta.json')))
        self.shards = self.meta['shards']
        self._db_cache = {}

    @staticmethod
    def _key_bytes(key_hash: int) -> bytes:
        return struct.pack('>Q', key_hash)

    def _shard_id(self, key_hash: int) -> int:
        return key_hash % self.shards

    def query(self, key_hash: int):
        import dbm
        sid = self._shard_id(key_hash)
        if sid not in self._db_cache:
            path = os.path.join(self.dir, f'seeds_{sid:03d}.dbm')
            self._db_cache[sid] = dbm.open(path, 'r')
        db = self._db_cache[sid]
        kb = self._key_bytes(key_hash)
        if kb in db:
            return pickle.loads(db[kb])  # List[(pid_i, pos)]
        return []

# ------------------------
# 位点倒排
# ------------------------
class LocusInverted:
    def __init__(self):
        self.index: Dict[str, dict] = {}
    @staticmethod
    def parse_locus_id(locus_id: str):
        m = re.match(r"([^:]+):([^:]+):(\d+)-(\d+):(.+)", locus_id)
        if not m: raise ValueError(f'bad locus_id: {locus_id}')
        return m.group(1), m.group(2), int(m.group(3)), int(m.group(4)), m.group(5)
    def build(self, gfa: GFA, pathseqs: Dict[str, PathSeq], anchor: int):
        pid_to_locus = dict(gfa.x_ptr)
        locus_to_pids = collections.defaultdict(list)
        for pid, lid in pid_to_locus.items():
            locus_to_pids[lid].append(pid)
        for locus_id, info in gfa.x_str.items():
            build, chrom, start, end, motif = self.parse_locus_id(locus_id)
            if chrom not in pathseqs: continue
            plen = len(pathseqs[chrom].seq)
            a0 = max(0, start-1-anchor)
            a1 = min(plen, end+anchor)
            self.index[locus_id] = {
                'motif': info['motif'],
                'candidate_Ls': info['candidate_Ls'],
                'chrom_path': chrom,
                'anchor': [a0, a1],
                'allele_paths': locus_to_pids.get(locus_id, []),
            }

# ------------------------
# 主流程：index
# ------------------------

def main():
    ap = argparse.ArgumentParser(description='PGG index（无SQL DBM 版）')
    ap.add_argument('--graph', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--method', choices=['syncmer','minimizer'], default='syncmer')
    ap.add_argument('--k', type=int, default=15)
    ap.add_argument('--s', type=int, default=5)
    ap.add_argument('--t', type=int, default=2)
    ap.add_argument('--w', type=int, default=10)
    ap.add_argument('--shards', type=int, default=64)
    ap.add_argument('--anchor', type=int, default=100)
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)

    # 读 GFA
    gfa = GFA(args.graph)

    # 路径扁平化
    pathseqs = PathFlattener.flatten(gfa)

    # 写路径序列池
    paths_dir = os.path.join(args.out, 'paths'); os.makedirs(paths_dir, exist_ok=True)
    for pid, ps in pathseqs.items():
        with open(os.path.join(paths_dir, f'{pid}.fa'), 'w') as w:
            w.write(f'>{pid}\n')
            s = ps.seq
            for i in range(0, len(s), 80):
                w.write(s[i:i+80] + '\n')

    # 路径坐标
    pickle.dump({pid: ps.segments for pid, ps in pathseqs.items()}, open(os.path.join(args.out, 'path_coord.pkl'), 'wb'), protocol=pickle.HIGHEST_PROTOCOL)

    # 位点倒排
    locus = LocusInverted(); locus.build(gfa, pathseqs, anchor=args.anchor)
    json.dump(locus.index, open(os.path.join(args.out, 'locus.json'), 'w'))

    # pid 映射（字符串 <-> 整数）
    pid_to_i = {pid:i for i, pid in enumerate(pathseqs.keys())}

    # 播种并写 DBM 分片
    seeder = Seeder()
    seeds_dir = os.path.join(args.out, 'seeds'); os.makedirs(seeds_dir, exist_ok=True)
    writer = DBMShardWriter(seeds_dir, shards=args.shards)
    writer.meta.update({'method': args.method, 'k': args.k, 's': args.s, 't': args.t, 'w': args.w})

    for pid, ps in pathseqs.items():
        pid_i = pid_to_i[pid]
        writer.add_path_meta(pid, len(ps.seq), pid_i)
        seq = ps.seq
        if not seq: continue
        if args.method == 'syncmer':
            poss = seeder.syncmer_positions(seq, args.k, args.s, args.t)
        else:
            poss = seeder.minimizer_positions(seq, args.k, args.w)
        for i in poss:
            kmer = Seeder.canonical(seq[i:i+args.k])
            kh = Seeder.kmer_hash(kmer)
            writer.add_seed(kh, pid_i, i)

    # 保存 pid 映射
    json.dump({'pid_to_i': pid_to_i}, open(os.path.join(args.out, 'pid_map.json'), 'w'))

    writer.finalize()
    print('[INDEX] done.')

if __name__ == '__main__':
    main()
