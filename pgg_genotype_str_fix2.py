#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
pgg_genotype_str.py  (JSONL 频率支持版，支持单端/双端 FASTQ)
---------------------------------------------------------
一个简化、可跑通的小型 STR 分型器（二倍体），支持：
- webSTR 风格 TSV 频率表（locus_id, L, AFR/EUR/.../ALL）
- JSONL/JSON 频率表（多种常见嵌套结构，见 load_freq_auto 注释）
- 混血先验 (--popmix "EUR=0.6,AFR=0.4")
- 极端先验 (--extreme "locus_id=L")
- 单端 FASTQ (--fq) 或双端 FASTQ (--fq1 + --fq2)

输出：每个位点一行 TSV，含 GT/GQ、like/post 两套分数、先验调试列、APL_json。
该脚本用于“验证先验是否进入打分链路”，非生产优化版本。
"""
import sys, os, gzip, json, math, argparse, csv, random
from collections import defaultdict

random.seed(1)

def _cap_gq_from_delta(delta):
    """把 Δ=post_best-post2 转成 GQ；对 inf/NaN 安全，并且封顶 99。"""
    if not math.isfinite(delta):
        return 99
    if delta < 0:
        delta = 0.0
    # 给个上限，避免数值过大（10 * 9.9 = 99）
    delta = min(delta, 9.9)
    return int(round(10.0 * delta))

# ------------------ FASTQ 读取 ------------------
def read_fastq(path):
    """读取单个 FASTQ(.gz)：yield (name, seq, qual)"""
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as f:
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

def load_reads_from_args(args):
    """
    根据命令行参数加载 reads：
      - 单端：--fq
      - 双端：--fq1 + --fq2

    双端模式下，R1 和 R2 中的每条 read 都视作独立 read 参与 SW 打分。
    为了区分，name 会在末尾补上 /1 或 /2（若原本没有的话）。
    """
    if args.fq and (args.fq1 or args.fq2):
        print("[ERROR] 不能同时使用 --fq 和 --fq1/--fq2，请二选一。", file=sys.stderr)
        sys.exit(1)
    if (args.fq1 and not args.fq2) or (args.fq2 and not args.fq1):
        print("[ERROR] 双端模式需要同时提供 --fq1 和 --fq2。", file=sys.stderr)
        sys.exit(1)
    if not args.fq and not (args.fq1 and args.fq2):
        print("[ERROR] 请使用 --fq（单端）或者 --fq1 + --fq2（双端）其一。", file=sys.stderr)
        sys.exit(1)

    reads = []
    if args.fq:
        # 单端
        print(f"[INFO] 使用单端 FASTQ：{args.fq}", file=sys.stderr)
        for name, seq, qual in read_fastq(args.fq):
            reads.append((name, seq, qual))
    else:
        # 双端
        print(f"[INFO] 使用双端 FASTQ：R1={args.fq1}  R2={args.fq2}", file=sys.stderr)
        # R1
        for name, seq, qual in read_fastq(args.fq1):
            if not name.endswith("/1"):
                name = name + "/1"
            reads.append((name, seq, qual))
        # R2
        for name, seq, qual in read_fastq(args.fq2):
            if not name.endswith("/2"):
                name = name + "/2"
            reads.append((name, seq, qual))

    if not reads:
        print("[ERROR] 未从 FASTQ 中读到任何 read，文件是否为空？", file=sys.stderr)
        sys.exit(2)
    print(f"[INFO] 共加载 reads 数量：{len(reads)}", file=sys.stderr)
    return reads

# ------------------ GFA 解析（X STR 行） ------------------
def parse_gfa_strs(path):
    loci = []  # list of dict
    with open(path, "r") as f:
        for ln in f:
            if not ln or ln[0] == "#":
                continue
            parts = ln.rstrip("\n").split("\t")
            if len(parts) >= 7 and parts[0] == "X" and parts[1] == "STR":
                locus_id, motif, left, right, LsCSV = parts[2], parts[3], parts[4], parts[5], parts[6]
                try:
                    cand_L = [int(x) for x in LsCSV.split(",") if x]
                except:
                    cand_L = []
                loci.append({
                    "locus_id": locus_id,
                    "motif": motif.upper(),
                    "left": left.upper(),
                    "right": right.upper(),
                    "candidates": cand_L,
                })
    return loci

# ------------------ 频率读取（自动检测 TSV / JSONL / JSON） ------------------
POPS = ["AFR","EUR","EAS","AMR","SAS","ALL"]

import csv, json

def load_freq_tsv(path):
    tab = defaultdict(dict)  # tab[locus_id][L] = {pop:val}
    with open(path, newline="", encoding="utf-8") as f:
        rd = csv.DictReader(f, delimiter="\t")
        for row in rd:
            lid = (row.get("locus_id") or row.get("id") or row.get("locus") or "").strip()
            if not lid:
                continue
            try:
                L = int(row.get("L"))
            except:
                continue
            d = {}
            for p in POPS:
                v = row.get(p, "")
                if v not in (None, "", "NA"):
                    try:
                        d[p] = float(v)
                    except:
                        pass
            if d:
                tab[lid][L] = d
    return tab

def _ingest_row_like(tab, obj):
    """尽可能鲁棒地从单个 JSON 对象 obj 提取频率。
    支持：
      A) 扁平: {"locus_id": "...", "L": 10, "ALL": 0.1, "EUR": 0.2, ...}
      B) pop->L: {"locus_id":"...", "freq": {"ALL":{"10":0.1,"11":0.2}, "EUR":{...}}}
      C) L->pop: {"locus_id":"...", "L_freq": {"10":{"ALL":0.1,"EUR":0.2}, "11":{...}}}
      D) 容器键别名: "priors"/"allele_freqs"/"pops"
    """
    lid = (obj.get("locus_id") or obj.get("pointer_id") or obj.get("id") or obj.get("locus") or "").strip()
    if not lid:
        return

    # A) 扁平
    if "L" in obj and any(k in obj for k in POPS):
        try:
            L = int(obj["L"])
            d = {}
            for p in POPS:
                if p in obj:
                    try:
                        d[p] = float(obj[p])
                    except:
                        pass
            if d:
                tab[lid][L] = d
                return
        except:
            pass

    # B/C/D) 容器
    cand_container_keys = ["freq","L_freq","priors","allele_freqs","pops"]
    container = None
    for k in cand_container_keys:
        if k in obj and isinstance(obj[k], dict):
            container = obj[k]
            break
    if container is None:
        # 也可能是 pop 键直接是 {L:val}
        has_pop_key = any((k in obj and isinstance(obj[k], dict)) for k in POPS)
        if has_pop_key:
            for p in POPS:
                if p in obj and isinstance(obj[p], dict):
                    for Ls, val in obj[p].items():
                        try:
                            L = int(Ls); val = float(val)
                            tab[lid].setdefault(L, {})[p] = val
                        except:
                            pass
            return
        # 或 obj 本身就是 {L:{pop:val}}
        try:
            for Ls, mp in obj.items():
                L = int(Ls)
                if isinstance(mp, dict):
                    for p, val in mp.items():
                        if p in POPS:
                            tab[lid].setdefault(L, {})[p] = float(val)
            return
        except:
            return
    else:
        # container 可能是 pop->L 或 L->pop
        if any(k in POPS for k in container.keys()):
            # pop->L
            for p, mp in container.items():
                if p not in POPS:
                    continue
                if isinstance(mp, dict):
                    it = mp.items()
                elif isinstance(mp, list):
                    it = []
                    for item in mp:
                        if isinstance(item, dict) and "L" in item:
                            it.append((item["L"], item.get("v", item.get("val", item.get("freq", None)))))
                else:
                    continue
                for Ls, val in it:
                    try:
                        L = int(Ls); val = float(val)
                        tab[lid].setdefault(L, {})[p] = val
                    except:
                        pass
        else:
            # L->pop
            for Ls, mp in container.items():
                try:
                    L = int(Ls)
                except:
                    continue
                if not isinstance(mp, dict):
                    continue
                for p, val in mp.items():
                    if p in POPS:
                        try:
                            tab[lid].setdefault(L, {})[p] = float(val)
                        except:
                            pass

def load_freq_json_lines(path, force_jsonl=False):
    tab = defaultdict(dict)

    if force_jsonl:
        # 严格 JSONL：逐行解析
        with open(path, "r", encoding="utf-8") as f:
            for ln in f:
                ln = ln.strip()
                if not ln or ln.startswith("#"):
                    continue
                try:
                    obj = json.loads(ln)
                    if isinstance(obj, (dict, list)):
                        if isinstance(obj, list):
                            for x in obj:
                                if isinstance(x, dict):
                                    _ingest_row_like(tab, x)
                        else:
                            _ingest_row_like(tab, obj)
                except:
                    continue
        return tab

    # AUTO：先尝试整体 JSON；失败则回退 JSONL
    try:
        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)   # 如果是 JSONL，这里多半会抛 Extra data
        # 成功：data 可能是 list 或 dict
        if isinstance(data, list):
            for x in data:
                if isinstance(x, dict):
                    _ingest_row_like(tab, x)
        elif isinstance(data, dict):
            _ingest_row_like(tab, data)
        return tab
    except Exception:
        # 回退 JSONL
        return load_freq_json_lines(path, force_jsonl=True)

def load_freq_auto(path):
    low = path.lower()
    if low.endswith(".jsonl"):
        return load_freq_json_lines(path, force_jsonl=True), "json"
    elif low.endswith(".json"):
        return load_freq_json_lines(path, force_jsonl=False), "json"
    else:
        return load_freq_tsv(path), "tsv"

# ------------------ 混血与极端先验 ------------------
def parse_popmix(s):
    # "EUR=0.6,AFR=0.4" -> dict
    if not s:
        return {}
    out = {}
    for kv in s.split(","):
        if not kv.strip():
            continue
        k,v = kv.split("=")
        out[k.strip().upper()] = float(v)
    # 归一化
    tot = sum(out.values())
    if tot > 0:
        for k in out:
            out[k] /= tot
    return out

def mix_prior_for_locus(freq_tab, locus_id, cand_L, popmix):
    # 返回 p_L（已平滑、归一化）
    eps = 1e-12
    raw = {}
    has_any = False
    for L in cand_L:
        pops = freq_tab.get(locus_id, {}).get(L, {})
        if popmix:
            s = 0.0
            for p,w in popmix.items():
                s += w * pops.get(p, 0.0)
            raw[L] = s
        else:
            # 默认 ALL
            raw[L] = pops.get("ALL", 0.0)
        has_any = has_any or (raw[L] > 0)
    # 平滑：若全为 0，则退为均匀
    if not has_any:
        for L in cand_L:
            raw[L] = 1.0
    # 加 eps 防零
    for L in cand_L:
        raw[L] = raw[L] + eps
    # 归一化
    tot = sum(raw.values())
    for L in cand_L:
        raw[L] /= tot
    return raw

def parse_extreme_list(lst):
    # ["locus_id=L", "locus2=7", ...]
    out = {}
    for s in lst or []:
        if "=" in s:
            lid, Ls = s.split("=",1)
            try:
                out[lid] = int(Ls)
            except:
                pass
    return out

def apply_extreme(prior_dict, extreme_L):
    # 把 prior 近似改成：L*=1-ε，其它=(ε/(n-1))
    if extreme_L is None:
        return prior_dict
    ks = list(prior_dict.keys())
    n = len(ks)
    if extreme_L not in ks or n<=1:
        return prior_dict
    eps = 1e-6
    for L in ks:
        prior_dict[L] = eps/(n-1)
    prior_dict[extreme_L] = 1.0 - eps
    return prior_dict

# ------------------ 简单双端 SW（线性gap；仅求最大分） ------------------
def sw_score(a, b, match=2, mismatch=-3, gap=-4, band=None):
    # 返回最大局部比对得分（不回溯）
    # 可选 band：若给定，限制 |i - j| <= band（节省时间）
    n, m = len(a), len(b)
    prev = [0]*(m+1)
    best = 0
    for i in range(1, n+1):
        cur = [0]*(m+1)
        ai = a[i-1]
        j_lo = 1
        j_hi = m
        if band is not None:
            j_lo = max(1, i - band)
            j_hi = min(m, i + band)
        for j in range(1, j_lo):
            cur[j] = 0
        for j in range(j_lo, j_hi+1):
            sc_diag = prev[j-1] + (match if ai == b[j-1] else mismatch)
            sc_up   = prev[j] + gap
            sc_left = cur[j-1] + gap
            sc = max(0, sc_diag, sc_up, sc_left)
            cur[j] = sc
            if sc > best:
                best = sc
        for j in range(j_hi+1, m+1):
            cur[j] = 0
        prev = cur
    return best

# ------------------ 基因型遍历与打分 ------------------
def all_genotype_pairs(cands):
    out = []
    for i, L1 in enumerate(cands):
        for j, L2 in enumerate(cands[i:]):
            out.append((L1, L2))
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gfa", required=True, help="包含 X\tSTR 行的 GFA 文件")
    # 单端或双端，二选一：
    ap.add_argument("--fq",  required=False, help="单端 FASTQ(.gz)")
    ap.add_argument("--fq1", required=False, help="双端 R1 FASTQ(.gz)")
    ap.add_argument("--fq2", required=False, help="双端 R2 FASTQ(.gz)")
    ap.add_argument("--freq", required=False, default=None, help="等位频率表，支持 .tsv / .jsonl / .json")
    ap.add_argument("--popmix", required=False, default=None, help='如 \"EUR=0.6,AFR=0.4\"，否则用 ALL 列或均匀')
    ap.add_argument("--extreme", required=False, default=None, nargs="*", help='形如 locus_id=L 的列表，用于极端先验验证')
    ap.add_argument("--tau", type=float, default=12.0, help="softmax 温度（把 SW 分数转成概率时的缩放）")
    ap.add_argument("--band", type=int, default=60, help="SW 带宽，限制 |i-j| 以提速")
    ap.add_argument("--max_reads_per_locus", type=int, default=500, help="每个位点最多使用多少条 read（随机下采样）")
    ap.add_argument("--out", required=True, help="输出 TSV 文件路径")
    args = ap.parse_args()

    loci = parse_gfa_strs(args.gfa)
    if not loci:
        print("GFA 中未找到 X STR 行", file=sys.stderr); sys.exit(2)

    freq_tab = {}
    freq_type = None
    if args.freq:
        freq_tab, freq_type = load_freq_auto(args.freq)
    popmix = parse_popmix(args.popmix)
    extreme_map = parse_extreme_list(args.extreme)

    # 读入全部 reads 到内存（小样本足够；双端即 R1+R2）
    reads = load_reads_from_args(args)

    out_rows = []
    for loc in loci:
        lid   = loc["locus_id"]
        motif = loc["motif"]
        left  = loc["left"]
        right = loc["right"]
        cands = sorted(set(loc["candidates"])) or list(range(3, 20))

        # H_L 序列池
        H = {L: (left + motif*L + right) for L in cands}

        # 抽样 reads
        idxs = list(range(len(reads)))
        random.shuffle(idxs)
        idxs = idxs[:args.max_reads_per_locus]
        reads_sel = [reads[i] for i in idxs]

        # 逐 read：把 SW 得分 softmax 成 P(r|L)
        per_read_log10PL = []
        for (qname, qseq, qqual) in reads_sel:
            # 计算每个 L 的 SW 得分
            scores = {}
            bestS = None
            for L, hseq in H.items():
                s = sw_score(qseq, hseq, match=2, mismatch=-3, gap=-4, band=args.band)
                scores[L] = s
                if (bestS is None) or (s > bestS):
                    bestS = s
            # softmax → 概率（以 log10 形式保存）
            log10P = {}
            denom = 0.0
            for L, s in scores.items():
                val = math.exp((s - bestS)/max(args.tau, 1e-6))
                log10P[L] = val
                denom += val
            for L in scores.keys():
                p = max(log10P[L]/denom, 1e-300)
                log10P[L] = math.log10(p)
            per_read_log10PL.append(log10P)

        n_reads = len(per_read_log10PL)

        # APL_json（按等位汇总的“代价”）：-10 * (sum_r log10 P(r|L) - max_L ...)
        sum_log10_P_L = {L: 0.0 for L in cands}
        for d in per_read_log10PL:
            for L,v in d.items():
                sum_log10_P_L[L] += v
        max_sum = max(sum_log10_P_L.values())
        APL = {str(L): int(round(-10.0*(sum_log10_P_L[L] - max_sum))) for L in cands}

        # 基因型穷举
        pairs = all_genotype_pairs(cands)

        # 先验（混血）；若有 extreme 指定，则覆盖
        prior = mix_prior_for_locus(freq_tab, lid, cands, popmix)
        if lid in extreme_map:
            prior = apply_extreme(prior, extreme_map[lid])

        like = {}
        post = {}
        prior_g = {}
        for (L1,L2) in pairs:
            # 似然
            s = 0.0
            for d in per_read_log10PL:
                pr = 0.5*(10**d[L1]) + 0.5*(10**d[L2])
                s += math.log10(max(pr, 1e-300))
            like[(L1,L2)] = s
            # 先验（HW）
            if L1 == L2:
                pg = prior[L1]*prior[L1]
            else:
                pg = 2.0*prior[L1]*prior[L2]
            prior_g[(L1,L2)] = pg
            post[(L1,L2)] = s + math.log10(max(pg, 1e-300))

        # 取最佳/次佳（like 与 post 各自独立排名）
        def top2(d):
            arr = sorted(d.items(), key=lambda kv: kv[1], reverse=True)
            best = arr[0]
            second = arr[1] if len(arr)>1 else (arr[0][0], float("-inf"))
            return best, second

        (best_like, sec_like) = top2(like)
        (best_post, sec_post) = top2(post)
        (gtL1, gtL2) = best_post[0]
        post_best, post2 = best_post[1], sec_post[1]
        like_best, like2 = best_like[1], sec_like[1]

        # 计算 Δ，若次优为 -inf，则视为 Δ=+inf → GQ=99
        delta = post_best - post2 if math.isfinite(post2) else float('inf')
        GQ = _cap_gq_from_delta(delta)

        # 计算 log_prior_best（基于最优基因型的 HW 先验）
        def _fmt(x):
            return x if math.isfinite(x) else "NA"

        row = {
            "locus_id": lid,
            "motif": motif,
            "n_reads": n_reads,
            "candidates": ",".join(str(x) for x in cands),
            "GT": f"{gtL1}/{gtL2}",
            "GQ": GQ,
            "like_best_log10": like_best,
            "like2_log10": _fmt(like2),
            "post_best_log10": post_best,
            "post2_log10": _fmt(post2),
            "log_prior_best": math.log10(max(prior_g[(gtL1,gtL2)], 1e-300)),
            "prior_src": (f"{freq_type or 'none'}|" + (("mix:"+args.popmix) if args.popmix else "ALL_or_uniform")),
            "prior_L_json": json.dumps({str(L): prior[L] for L in cands}),
            "APL_json": json.dumps(APL, ensure_ascii=False),
        }
        out_rows.append(row)

    # 写出 TSV
    cols = ["locus_id","motif","n_reads","candidates","GT","GQ",
            "like_best_log10","like2_log10","post_best_log10","post2_log10",
            "log_prior_best","prior_src","prior_L_json","APL_json"]
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in out_rows:
            w.writerow(r)

if __name__ == "__main__":
    main()
