from __future__ import annotations
import os
import sys
import gzip
import json
import struct
import argparse
import collections
import hashlib
import multiprocessing
from dataclasses import dataclass
from typing import Iterable, Tuple, List, Dict
import pickle
from functools import partial

# Try importing tqdm, otherwise use a simple alternative
try:
    from tqdm import tqdm
except ImportError:
    def tqdm(iterable, desc="", total=None):
        for item in iterable:
            yield item

# ------------------------
# General Utilities
# ------------------------
DNA_COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')

def revcomp(seq: str) -> str:
    return seq.translate(DNA_COMP)[::-1]

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
        if k > n or s > k or t < 0 or t > k - s:
            return []
        pos = []
        for i in range(0, n - k + 1):
            kmer = seq[i:i + k]
            minv = None
            minp = None
            for j in range(0, k - s + 1):
                v = kmer[j:j + s]
                if (minv is None) or (v < minv):
                    minv, minp = v, j
            if minp == t:
                pos.append(i)
        return pos

    @staticmethod
    def minimizer_positions(seq: str, k: int, w: int):
        n = len(seq)
        if k > n or w <= 0:
            return []
        pos = []
        window = []  # (hash, i)

        def kmer(i):
            rc = revcomp(seq[i:i + k])
            s = seq[i:i + k]
            return s if s <= rc else rc

        for i in range(0, n - k + 1):
            import hashlib as _h
            h = int.from_bytes(
                _h.blake2b(kmer(i).encode('ascii'), digest_size=8).digest(),
                'big'
            )
            window.append((h, i))
            left = i - w + 1
            while window and window[0][1] < left:
                window.pop(0)
            if left >= 0:
                m = min(window, key=lambda x: x[0])
                if not pos or pos[-1] != m[1]:
                    pos.append(m[1])
        return pos

# ------------------------
# Index Reading (DBM Shards)
# ------------------------
class DBMShardReader:
    def __init__(self, dir_path: str):
        import dbm
        self.dir = dir_path
        self.meta = json.load(open(os.path.join(self.dir, 'meta.json')))
        self.shards = self.meta['shards']
        self._db_cache = {}
        self.pid_to_i: Dict[str, int] = self.meta.get('pid_to_i', {})
        self.i_to_pid: Dict[int, str] = {
            int(k): v for k, v in self.meta.get('i_to_pid', {}).items()
        }
        self.path_len: Dict[str, int] = self.meta.get('paths', {})

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
            # Open in 'r' read-only mode, usually supports concurrent reading
            try:
                self._db_cache[sid] = dbm.open(path, 'r')
            except Exception as e:
                # Compatibility: some dbm implementations might need 'c' or other flags,
                # but read-only should use 'r'.
                # If multiprocessing, must ensure opening within the process.
                raise RuntimeError(f"Cannot open DBM file {path}: {e}")
                
        db = self._db_cache[sid]
        kb = self._key_bytes(key_hash)
        if kb in db:
            return pickle.loads(db[kb])  # List[(pid_i, pos)]
        return []

# ------------------------
# FASTQ Reading: Iterator mode (removed pre-counting of lines)
# ------------------------

def read_fastq(path):
    """Read single FASTQ(.gz): yield (name, seq, qual)"""
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt", encoding='utf-8', errors='replace') as f:
        while True:
            name = f.readline().rstrip()
            if not name:
                break
            seq  = f.readline().rstrip()
            plus = f.readline().rstrip()
            qual = f.readline().rstrip()
            if not qual:
                break
            yield name[1:] if name.startswith("@") else name, seq, qual

def chunked_iterable(iterable, size):
    """Split iterator into chunks for multiprocessing distribution"""
    import itertools
    it = iter(iterable)
    while True:
        chunk = tuple(itertools.islice(it, size))
        if not chunk:
            break
        yield chunk

# ------------------------
# Chaining and GAF Output
# ------------------------

@dataclass
class Hit:
    pid: str
    ppos: int
    rpos: int

class Mapper:
    def __init__(self, idx_dir: str):
        self.seeds = DBMShardReader(os.path.join(idx_dir, 'seeds'))
        self.meta = self.seeds.meta
        self.path_len = self.meta['paths']
        pid_map = json.load(open(os.path.join(idx_dir, 'pid_map.json')))
        self.i_to_pid = {i: pid for pid, i in pid_map['pid_to_i'].items()}

    def seed_hits(
        self,
        seq: str,
        method: str,
        k: int,
        s: int,
        t: int,
        w: int,
        max_occ: int = 200
    ) -> List[Hit]:
        hits: List[Hit] = []
        strands = [(seq, '+'), (revcomp(seq), '-')]
        for qseq, strand in strands:
            if method == 'syncmer':
                qpos = Seeder.syncmer_positions(qseq, k, s, t)
            else:
                qpos = Seeder.minimizer_positions(qseq, k, w)
            for qp in qpos:
                kmer = Seeder.canonical(qseq[qp:qp + k])
                kh = Seeder.kmer_hash(kmer)
                postings = self.seeds.query(kh)
                if not postings or (len(postings) > max_occ):
                    continue
                for pid_i, ppos in postings:
                    ridx = qp if strand == '+' else (len(seq) - k - qp)
                    hits.append(
                        Hit(pid=self.i_to_pid[pid_i], ppos=ppos, rpos=ridx)
                    )
        return hits

    @staticmethod
    def chain(hits: List[Hit], slack: int = 10):
        by_path = collections.defaultdict(list)
        for h in hits:
            by_path[h.pid].append(h)
        best_pid, best_chain = None, []
        for pid, hs in by_path.items():
            buckets = collections.defaultdict(list)
            for h in hs:
                d = h.ppos - h.rpos
                key = int(round(d / slack))
                buckets[key].append(h)
            key, chain_hits = max(
                buckets.items(), key=lambda kv: len(kv[1])
            )
            chain_hits.sort(key=lambda x: x.ppos)
            if len(chain_hits) > len(best_chain):
                best_pid, best_chain = pid, chain_hits
        return best_pid, best_chain

    def gaf_line(self, rid: str, rseq: str, pid: str, chain: List[Hit], k: int) -> str:
        if not chain or pid is None:
            return f"{rid}\t{len(rseq)}\t0\t0\t+\t*\t0\t0\t0\t0\t0\tcg:Z:*"
        qstart = min(h.rpos for h in chain)
        qend   = max(h.rpos for h in chain) + k
        pstart = min(h.ppos for h in chain)
        pend   = max(h.ppos for h in chain) + k
        plen   = self.path_len.get(pid, 0)
        nm = 0
        mapq = min(60, 10 + 5 * len(chain))
        cg = f"SEEDS:{len(chain)}"
        return (
            f"{rid}\t{len(rseq)}\t{qstart}\t{qend}\t+\t{pid}\t{plen}\t"
            f"{pstart}\t{pend}\t{nm}\t{mapq}\tcg:Z:{cg}"
        )

# ------------------------
# Multiprocessing Worker Logic
# ------------------------

# Global variable for the worker process to hold the Mapper instance, avoiding repeated loading
worker_mapper = None

def init_worker(index_dir):
    """Worker process initialization: load index only once per process"""
    global worker_mapper
    worker_mapper = Mapper(index_dir)

def process_batch_single(batch, args_dict):
    """Process a batch of single-end reads"""
    global worker_mapper
    results = []
    method, k, s, t, w, max_occ, slack = args_dict['method'], args_dict['k'], args_dict['s'], args_dict['t'], args_dict['w'], args_dict['max_occ'], args_dict['slack']
    
    for rid, seq, _ in batch:
        hits = worker_mapper.seed_hits(seq, method, k, s, t, w, max_occ)
        pid, chain = worker_mapper.chain(hits, slack)
        line = worker_mapper.gaf_line(rid, seq, pid, chain, k)
        results.append(line)
    return results

def process_batch_paired(batch, args_dict):
    """Process a batch of paired-end reads"""
    global worker_mapper
    results = []
    method, k, s, t, w, max_occ, slack = args_dict['method'], args_dict['k'], args_dict['s'], args_dict['t'], args_dict['w'], args_dict['max_occ'], args_dict['slack']
    
    for (h1, s1, _), (h2, s2, _) in batch:
        # Process R1
        hits1 = worker_mapper.seed_hits(s1, method, k, s, t, w, max_occ)
        pid1, chain1 = worker_mapper.chain(hits1, slack)
        rid1_out = h1 if h1.endswith('/1') else h1 + '/1'
        line1 = worker_mapper.gaf_line(rid1_out, s1, pid1, chain1, k)
        
        # Process R2
        hits2 = worker_mapper.seed_hits(s2, method, k, s, t, w, max_occ)
        pid2, chain2 = worker_mapper.chain(hits2, slack)
        rid2_out = h2 if h2.endswith('/2') else h2 + '/2'
        line2 = worker_mapper.gaf_line(rid2_out, s2, pid2, chain2, k)
        
        results.append(line1 + "\n" + line2)
    return results

# ------------------------
# CLI
# ------------------------

def main():
    ap = argparse.ArgumentParser(description='PGG map (Multiprocessing Optimized)')
    ap.add_argument('--index', required=True, help='Index directory')
    ap.add_argument('--reads', help='Single-end FASTQ(.gz)')
    ap.add_argument('--reads1', help='Paired-end R1 FASTQ(.gz)')
    ap.add_argument('--reads2', help='Paired-end R2 FASTQ(.gz)')
    ap.add_argument('--out', required=True, help='Output GAF file')
    
    # Algorithm parameters
    ap.add_argument('--method', choices=['syncmer', 'minimizer'], default='syncmer')
    ap.add_argument('--k', type=int, default=15)
    ap.add_argument('--s', type=int, default=5)
    ap.add_argument('--t', type=int, default=2, help="Syncmer parameter (Do NOT enter CPU thread count!)")
    ap.add_argument('--w', type=int, default=10)
    ap.add_argument('--max_occ', type=int, default=200)
    ap.add_argument('--slack', type=int, default=10)
    
    # New performance parameters
    ap.add_argument('--threads', type=int, default=1, help='Number of CPU threads for parallel processing')
    ap.add_argument('--batch_size', type=int, default=5000, help='Number of reads processed per thread at a time')

    args = ap.parse_args()

    # Pack parameters
    algo_params = {
        'method': args.method, 'k': args.k, 's': args.s, 't': args.t, 
        'w': args.w, 'max_occ': args.max_occ, 'slack': args.slack
    }

    # Prepare input streams
    if args.reads:
        mode = 'single'
        print(f"[MAP] Single-end mode, input: {args.reads}", file=sys.stderr)
        read_iter = read_fastq(args.reads)
    elif args.reads1 and args.reads2:
        mode = 'paired'
        print(f"[MAP] Paired-end mode, input: R1={args.reads1} R2={args.reads2}", file=sys.stderr)
        # Use zip to read both files simultaneously
        read_iter = zip(read_fastq(args.reads1), read_fastq(args.reads2))
    else:
        print("[ERROR] Please specify --reads or (--reads1 + --reads2)", file=sys.stderr)
        sys.exit(1)

    print(f"[MAP] Starting multiprocessing: {args.threads} threads (Batch size: {args.batch_size})", file=sys.stderr)

    # Start multiprocessing pool
    # initializer ensures each child process opens DBM index independently to avoid lock conflicts
    with multiprocessing.Pool(processes=args.threads, initializer=init_worker, initargs=(args.index,)) as pool:
        
        # Split reads generator into chunks
        chunks = chunked_iterable(read_iter, args.batch_size)
        
        # Select processing function
        worker_func = process_batch_single if mode == 'single' else process_batch_paired
        func_with_args = partial(worker_func, args_dict=algo_params)

        with open(args.out, 'w') as outf:
            # imap provides ordered results and saves some memory
            # Use tqdm to show the number of processed batches
            for batch_results in tqdm(pool.imap(func_with_args, chunks), desc="Processing Batches"):
                for line in batch_results:
                    outf.write(line + "\n")

    print('[MAP] done.', file=sys.stderr)

if __name__ == '__main__':
    main()
