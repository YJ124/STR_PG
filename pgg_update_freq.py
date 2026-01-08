#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
pgg_update_freq.py

根据单个样本的 STR 分型结果 (genotypes.tsv)，
对种群频率库 freq.jsonl 做“增量更新”（只使用纯合位点）。

逻辑：
  - 只统计 GT 为 x/x 的纯合基因型；
  - 每个纯合样本对该等位基因贡献 2 个拷贝；
  - 在 freq.jsonl 中为对应 pointer_id + 人群 标签更新：
        counts[pop][allele] += 2
        total_alleles[pop]  += 2
    然后重新计算：
        freq[pop][allele] = counts[pop][allele] / total_alleles[pop]

说明：
  - 为兼容已有 freq.jsonl，脚本会：
      * 若原记录没有 counts/total_alleles 字段，则用 freq * init_total_alleles
        推出一个“伪计数”，再继续增量更新。
  - locus_id 与 pointer_id 默认认为是同一个字符串（与你现有 pipeline 一致）。

使用示例：

  python pgg_update_freq.py \\
      --geno genotypes.tsv \\
      --freq-in freq19.jsonl \\
      --freq-out freq19.updated.jsonl \\
      --pop EAS \\
      --min-GQ 20 \\
      --init-total-alleles 1000
"""

import os
import json
import argparse
from typing import Dict, Any, Tuple


def load_freq_jsonl(path: str,
                    init_total_alleles: int) -> Dict[str, Dict[str, Any]]:
    """
    读取 freq.jsonl，返回：
        { pointer_id: { 'pointer_id':..., 'freq':..., 'counts':..., 'total_alleles':... } }
    若原始 jsonl 中没有 counts/total_alleles，则用 freq * init_total_alleles 生成伪计数。
    """
    db: Dict[str, Dict[str, Any]] = {}
    if not os.path.exists(path):
        return db

    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            obj = json.loads(line)
            pid = obj.get("pointer_id")
            if not pid:
                continue

            freq = obj.get("freq", {})
            counts = obj.get("counts", {})
            total_alleles = obj.get("total_alleles", {})

            # 如果没有 counts/total_alleles，用 freq 推一个初始计数出来
            if not counts or not total_alleles:
                counts = {}
                total_alleles = {}
                for pop, dist in freq.items():
                    pop_counts = {}
                    for L_str, p in dist.items():
                        # 用频率 * init_total_alleles 得到一个近似整数计数
                        c = max(0, int(round(float(p) * init_total_alleles)))
                        if c > 0:
                            pop_counts[L_str] = c
                    total_c = sum(pop_counts.values())
                    # 如果某个人群没有任何等位计数，则 total_alleles 为 0
                    counts[pop] = pop_counts
                    total_alleles[pop] = total_c

            obj["freq"] = freq
            obj["counts"] = counts
            obj["total_alleles"] = total_alleles
            db[pid] = obj

    return db


def save_freq_jsonl(db: Dict[str, Dict[str, Any]], out_path: str) -> None:
    """
    将更新后的频率库写回 JSONL 文件。保持 pointer_id/freq/counts/total_alleles 四个字段。
    """
    with open(out_path, "w", encoding="utf-8") as f:
        for pid, obj in db.items():
            # 确保 pointer_id 一致
            obj["pointer_id"] = pid
            # 先根据 counts & total_alleles 重新计算 freq
            freq = obj.get("freq", {})
            counts = obj.get("counts", {})
            total_alleles = obj.get("total_alleles", {})

            new_freq = {}
            for pop, alle_counts in counts.items():
                ta = total_alleles.get(pop, 0)
                if ta <= 0:
                    continue
                pop_freq = {}
                for L_str, c in alle_counts.items():
                    if c <= 0:
                        continue
                    pop_freq[L_str] = c / float(ta)
                if pop_freq:
                    new_freq[pop] = pop_freq
            obj["freq"] = new_freq

            f.write(json.dumps(obj, ensure_ascii=False) + "\n")


def parse_gt(gt: str) -> Tuple[int, int]:
    """
    将 GT 形如 '10/12' 解析成 (10, 12)。
    若为 './.' 或不合法，抛出异常由上层处理。
    """
    if gt is None:
        raise ValueError("GT is None")
    gt = gt.strip()
    if gt == "./.":
        raise ValueError("missing GT")
    if "/" not in gt:
        raise ValueError(f"bad GT format: {gt}")
    a, b = gt.split("/", 1)
    return int(a), int(b)


def update_freq_with_sample(
    db: Dict[str, Dict[str, Any]],
    geno_tsv: str,
    pop: str,
    min_gq: int = 0
) -> None:
    """
    核心：用一个样本的 genotypes.tsv 更新频率库 db。
    只计入：
      - GT 纯合（L/L）
      - GQ >= min_gq
    """
    if pop is None or pop == "":
        raise ValueError("pop（人群标签）不能为空，例如 EAS/EUR/AFR 等。")

    with open(geno_tsv, "r", encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split("\t")
        # 期待的列名（与你原始 pgg_genotype.tsv 保持一致）：
        # locus_id, motif, n_reads, candidates, GT, GQ, post_best_log10, post2_log10, APL_json
        col_idx = {name: i for i, name in enumerate(header)}

        required_cols = ["locus_id", "GT", "GQ"]
        for col in required_cols:
            if col not in col_idx:
                raise RuntimeError(
                    f"TSV 中缺少必须列: {col}，当前列: {header}"
                )

        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            locus_id = parts[col_idx["locus_id"]]
            gt_str = parts[col_idx["GT"]]
            gq_str = parts[col_idx["GQ"]]

            # 跳过缺失 GT
            if gt_str == "./.":
                continue

            # GQ 过滤
            try:
                gq = int(gq_str)
            except Exception:
                gq = 0
            if gq < min_gq:
                continue

            # 解析 GT
            try:
                L1, L2 = parse_gt(gt_str)
            except Exception:
                continue

            # 只计入纯合位点
            if L1 != L2:
                continue
            L = L1

            # 开始对频率库做增量更新
            pid = locus_id  # 默认：locus_id == pointer_id
            if pid not in db:
                # 该位点在原 freq.jsonl 中不存在，新建一条记录
                db[pid] = {
                    "pointer_id": pid,
                    "freq": {},
                    "counts": {},
                    "total_alleles": {}
                }

            obj = db[pid]
            freq = obj.setdefault("freq", {})
            counts = obj.setdefault("counts", {})
            total_alleles = obj.setdefault("total_alleles", {})

            pop_counts = counts.setdefault(pop, {})
            pop_total = total_alleles.get(pop, 0)

            L_str = str(L)

            # 对纯合个体：贡献 2 个该等位基因拷贝
            c_old = pop_counts.get(L_str, 0)
            c_new = c_old + 2
            pop_counts[L_str] = c_new
            pop_total += 2
            total_alleles[pop] = pop_total

            # freq 先暂时不算，这里只修改 counts / total_alleles，
            # 最终写出时再统一用 counts/total_alleles 计算 freq。
            obj["counts"] = counts
            obj["total_alleles"] = total_alleles
            db[pid] = obj


def main():
    ap = argparse.ArgumentParser(
        description="根据 genotypes.tsv 的纯合位点对 freq.jsonl 做增量更新"
    )
    ap.add_argument("--geno", required=True,
                    help="pgg_genotype 输出的 genotypes.tsv 路径")
    ap.add_argument("--freq-in", required=True,
                    help="原始 freq.jsonl / freq19.jsonl 路径")
    ap.add_argument("--freq-out", required=True,
                    help="更新后的 freq.jsonl 输出路径")
    ap.add_argument("--pop", required=True,
                    help="该样本所属人群标签（如 EAS/EUR/AFR 等）")
    ap.add_argument("--min-GQ", type=int, default=0,
                    help="只使用 GQ >= 此阈值 的位点（默认 0）")
    ap.add_argument("--init-total-alleles", type=int, default=1000,
                    help="当原始 freq.jsonl 没有 counts/total_alleles 时，"
                         "用 freq * init_total_alleles 生成伪计数的基数（默认 1000）")

    args = ap.parse_args()

    # 1) 读取旧频率库
    print(f"[INFO] loading freq jsonl from: {args.freq_in}")
    db = load_freq_jsonl(args.freq_in, init_total_alleles=args.init_total_alleles)
    print(f"[INFO] loaded {len(db)} loci from freq jsonl.")

    # 2) 用该样本的 genotypes.tsv 做增量更新
    print(f"[INFO] updating freq with genotypes from: {args.geno}")
    update_freq_with_sample(
        db,
        geno_tsv=args.geno,
        pop=args.pop,
        min_gq=args.min_GQ
    )

    # 3) 写出新的 freq.jsonl
    print(f"[INFO] writing updated freq jsonl to: {args.freq_out}")
    save_freq_jsonl(db, args.freq_out)
    print("[DONE] freq jsonl updated.")


if __name__ == "__main__":
    main()
