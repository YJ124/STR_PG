from __future__ import annotations

import csv
import json
import math
from dataclasses import asdict
from itertools import combinations_with_replacement
from pathlib import Path
from typing import Iterable

import numpy as np

from .alignment import AlignmentSoftmaxParameters, alignment_softmax_log_likelihood_matrix
from .defaults import DEFAULT_LIKELIHOOD_BACKEND
from .io import iter_fastq, read_gaf
from .likelihoods import create_likelihood_backend
from .mixture import fit_allele_mixture
from .models import AlleleRecord, GenotypeCall, LocusRecord, ProvisionalAllele
from .novel import candidate_repeat_counts, has_anchored_repeat_evidence
from .phmm import MotifAwarePHMM, PHMMParameters
from .priors import genotype_priors, parse_population_mix
from .registry import AlleleRegistry
from .read_selection import (
    select_informative_read_records,
    selected_reads_by_locus,
    write_read_cache,
)
from .utils import log_normalize, logsumexp, open_text, phred_from_error_probability


def collect_locus_reads(
    gaf_path: str | Path,
    reads1: str | Path | None = None,
    reads2: str | Path | None = None,
    min_mapq: int = 0,
    bam_path: str | Path | None = None,
) -> dict[str, list[str]]:
    assignments: dict[tuple[str, int], str] = {}
    for aln in read_gaf(gaf_path):
        if aln.mapq < min_mapq:
            continue
        assignments[(aln.query_name, aln.mate)] = aln.pointer_id
    by_locus: dict[str, list[str]] = {}
    if bam_path is not None:
        from .io import iter_bam_reads
        for name, seq, mate in iter_bam_reads(bam_path):
            pid = assignments.get((name, mate)) or assignments.get((name, 0))
            if pid:
                by_locus.setdefault(pid, []).append(seq)
        return by_locus
    if reads1 is None:
        raise ValueError("Provide FASTQ reads1 or BAM/CRAM input")
    for name, seq, _ in iter_fastq(reads1):
        pid = assignments.get((name, 1)) or assignments.get((name, 0))
        if pid:
            by_locus.setdefault(pid, []).append(seq)
    if reads2:
        for name, seq, _ in iter_fastq(reads2):
            pid = assignments.get((name, 2))
            if pid:
                by_locus.setdefault(pid, []).append(seq)
    return by_locus


def allele_log_likelihood_matrix(
    reads: list[str],
    alleles: list[AlleleRecord],
    phmm: MotifAwarePHMM,
    likelihood_model: str = DEFAULT_LIKELIHOOD_BACKEND,
    alignment_params: AlignmentSoftmaxParameters | None = None,
) -> np.ndarray:
    normalized = "sw" if likelihood_model == "alignment_softmax" else likelihood_model
    backend = create_likelihood_backend(
        normalized,
        sw_params=alignment_params,
        phmm_params=phmm.params,
    )
    return backend.score_reads(reads, alleles)


def _score_genotypes(
    log_likelihoods: np.ndarray,
    alleles: list[AlleleRecord],
    priors: dict[tuple[int, int], float],
) -> list[dict]:
    results: list[dict] = []
    for i, j in combinations_with_replacement(range(len(alleles)), 2):
        per_read = np.array([logsumexp([row[i], row[j]]) - math.log(2.0) for row in log_likelihoods])
        log_likelihood = float(per_read.sum())
        log_prior = math.log(max(priors[(i, j)], 1e-300))
        results.append({"i": i, "j": j, "log_likelihood": log_likelihood, "log_prior": log_prior, "log_score": log_likelihood + log_prior})
    results.sort(key=lambda row: row["log_score"], reverse=True)
    posteriors = log_normalize([row["log_score"] for row in results])
    for row, posterior in zip(results, posteriors):
        row["posterior"] = posterior
    return results


def call_locus(
    locus: LocusRecord,
    reads: list[str],
    phmm: MotifAwarePHMM,
    population_mix: dict[str, float],
    theta: float,
    sigma: float,
    heterogeneity_threshold: float | None = None,
    alleles_override: list[AlleleRecord] | None = None,
    prior_mode: str = "full",
    use_length_smoothing: bool = True,
    likelihood_model: str = DEFAULT_LIKELIHOOD_BACKEND,
    alignment_params: AlignmentSoftmaxParameters | None = None,
    use_mixture_model: bool = True,
    precomputed_log_likelihoods: np.ndarray | None = None,
) -> tuple[GenotypeCall, np.ndarray, list[dict]]:
    alleles = alleles_override or locus.alleles
    if not reads:
        raise ValueError("Cannot genotype locus without reads")
    if precomputed_log_likelihoods is None:
        matrix = allele_log_likelihood_matrix(
            reads, alleles, phmm, likelihood_model, alignment_params
        )
    else:
        matrix = np.asarray(precomputed_log_likelihoods, dtype=float)
        if matrix.shape != (len(reads), len(alleles)):
            raise ValueError(
                "Precomputed likelihood matrix shape does not match reads/alleles: "
                f"{matrix.shape} != {(len(reads), len(alleles))}"
            )
    priors = genotype_priors(
        alleles,
        population_mix,
        theta=theta,
        sigma=sigma,
        prior_mode=prior_mode,
        use_length_smoothing=use_length_smoothing,
    )
    genotype_rows = _score_genotypes(matrix, alleles, priors)
    best = genotype_rows[0]
    posterior = float(best["posterior"])
    gq = phred_from_error_probability(1.0 - posterior)

    if use_mixture_model:
        mixture = fit_allele_mixture(matrix)
        best_likelihood_only = max(row["log_likelihood"] for row in genotype_rows)
        log_mix: float | None = mixture.log_likelihood
        d_paper: float | None = max(
            0.0, -2.0 * (best["log_score"] - mixture.log_likelihood)
        )
        d_likelihood: float | None = max(
            0.0, -2.0 * (best_likelihood_only - mixture.log_likelihood)
        )
        weights: dict[int, float] | None = {
            alleles[idx].repeat_count: float(weight)
            for idx, weight in enumerate(mixture.weights)
            if weight >= 1e-4
        }
    else:
        log_mix = None
        d_paper = None
        d_likelihood = None
        weights = None
    a1, a2 = alleles[best["i"]], alleles[best["j"]]
    if heterogeneity_threshold is None:
        flag = "." if use_mixture_model else "NA"
    elif not use_mixture_model:
        flag = "NA"
    else:
        assert d_paper is not None
        flag = "1" if d_paper >= heterogeneity_threshold else "0"
    call = GenotypeCall(
        pointer_id=locus.pointer_id,
        chrom=locus.chrom,
        start=locus.start,
        end=locus.end,
        motif=locus.motif,
        depth=len(reads),
        allele1=a1.repeat_count,
        allele2=a2.repeat_count,
        gt=f"{a1.repeat_count}/{a2.repeat_count}",
        gq=gq,
        posterior=posterior,
        log_likelihood=float(best["log_likelihood"]),
        log_prior=float(best["log_prior"]),
        log_score=float(best["log_score"]),
        log_mix=log_mix,
        d_paper=d_paper,
        d_likelihood=d_likelihood,
        mixture_weights=weights,
        heterogeneity_threshold=heterogeneity_threshold,
        heterogeneity_flag=flag,
    )
    return call, matrix, genotype_rows


def _temporary_novel_allele(locus: LocusRecord, repeat_count: int) -> AlleleRecord:
    sequence = locus.left_flank + locus.motif * repeat_count + locus.right_flank
    populations = {pop for allele in locus.alleles for pop in allele.frequencies}
    return AlleleRecord(
        allele_id=f"{locus.pointer_id}:NOVEL:L{repeat_count}",
        repeat_count=repeat_count,
        sequence=sequence,
        repeat_sequence=locus.motif * repeat_count,
        motif=locus.motif,
        repeat_start=len(locus.left_flank),
        repeat_end=len(locus.left_flank) + len(locus.motif) * repeat_count,
        frequencies={pop: 1e-6 for pop in populations or {"GLOBAL"}},
        counts={pop: 0.0 for pop in populations or {"GLOBAL"}},
        status="provisional",
        source="read_hypothesis",
    )


def discover_novel_alleles(
    locus: LocusRecord,
    reads: list[str],
    phmm: MotifAwarePHMM,
    population_mix: dict[str, float],
    theta: float,
    sigma: float,
    baseline_call: GenotypeCall,
    min_support: int = 3,
    min_log_likelihood_gain: float = 3.0,
    min_gq: int = 20,
    anchor_k: int = 12,
    stutter_max_fraction: float = 0.15,
    likelihood_backend: str = DEFAULT_LIKELIHOOD_BACKEND,
    alignment_params: AlignmentSoftmaxParameters | None = None,
    prior_mode: str = "full",
    use_length_smoothing: bool = True,
) -> tuple[list[ProvisionalAllele], list[AlleleRecord]]:
    hypotheses: list[ProvisionalAllele] = []
    accepted_alleles = list(locus.alleles)
    candidates = candidate_repeat_counts(reads, locus, min_support=min_support, anchor_k=anchor_k)
    registered_counts = [a.repeat_count for a in locus.alleles]
    for repeat_count, support in sorted(candidates.items(), key=lambda kv: (-kv[1], kv[0])):
        nearest_distance = min(abs(repeat_count - x) for x in registered_counts) if registered_counts else 999
        stutter_shadow = nearest_distance == 1 and support / max(1, len(reads)) <= stutter_max_fraction
        novel = _temporary_novel_allele(locus, repeat_count)
        extended = accepted_alleles + [novel]
        call, _, _ = call_locus(
            locus,
            reads,
            phmm,
            population_mix,
            theta,
            sigma,
            alleles_override=extended,
            likelihood_model=likelihood_backend,
            alignment_params=alignment_params,
            prior_mode=prior_mode,
            use_length_smoothing=use_length_smoothing,
        )
        gain = call.log_score - baseline_call.log_score
        reasons: list[str] = []
        if stutter_shadow:
            reasons.append("one-step stutter shadow")
        if gain < min_log_likelihood_gain:
            reasons.append(f"log-score gain {gain:.3f} < {min_log_likelihood_gain}")
        if call.gq < min_gq:
            reasons.append(f"GQ {call.gq} < {min_gq}")
        if support < min_support:
            reasons.append(f"support {support} < {min_support}")
        accepted = not reasons and repeat_count in {call.allele1, call.allele2}
        hypothesis = ProvisionalAllele(
            repeat_count=repeat_count,
            sequence=novel.sequence,
            support=support,
            flank_support=support,
            log_likelihood_gain=gain,
            gq=call.gq,
            status="accepted" if accepted else "provisional",
            reason="all acceptance criteria passed" if accepted else "; ".join(reasons or ["novel allele not selected in MAP genotype"]),
        )
        hypotheses.append(hypothesis)
        if accepted:
            accepted_alleles.append(novel)
            baseline_call = call
    return hypotheses, accepted_alleles


OUTPUT_FIELDS = [
    "pointer_id", "chrom", "start", "end", "motif", "mapped_depth", "depth", "GT", "GQ", "allele1", "allele2",
    "posterior", "log_likelihood", "log_prior", "log_score", "log_mix", "D_PAPER", "D_LIKELIHOOD",
    "heterogeneity_threshold", "HET_FLAG", "mixture_weights",
]


def genotype_dataset(
    registry_path: str | Path,
    gaf_path: str | Path,
    reads1: str | Path | None,
    out_path: str | Path,
    reads2: str | Path | None = None,
    population_mix_text: str | None = None,
    theta: float = 0.01,
    sigma: float = 10.0,
    min_mapq: int = 0,
    phmm_params: PHMMParameters | None = None,
    heterogeneity_threshold: float | None = None,
    discover_novel: bool = False,
    novel_out: str | Path | None = None,
    min_novel_support: int = 3,
    min_novel_log_gain: float = 3.0,
    min_novel_gq: int = 20,
    bam_path: str | Path | None = None,
    require_anchored_repeat: bool = True,
    min_reference_log_gain: float | None = None,
    prior_mode: str = "full",
    use_length_smoothing: bool = True,
    likelihood_backend: str = DEFAULT_LIKELIHOOD_BACKEND,
    likelihood_model: str | None = None,
    alignment_params: AlignmentSoftmaxParameters | None = None,
    use_mixture_model: bool = True,
    sample_id: str = "sample",
    read_cache_dir: str | Path | None = None,
    use_read_cache: bool = True,
) -> dict[str, int | str]:
    registry = AlleleRegistry.load(registry_path)
    population_mix = parse_population_mix(population_mix_text)
    if likelihood_model is not None:
        alias = "sw" if likelihood_model == "alignment_softmax" else likelihood_model
        if likelihood_backend != DEFAULT_LIKELIHOOD_BACKEND and alias != likelihood_backend:
            raise ValueError(
                "Conflicting likelihood_backend and deprecated likelihood_model"
            )
        likelihood_backend = alias
    if likelihood_backend not in {"sw", "phmm"}:
        raise ValueError("likelihood_backend must be 'sw' or 'phmm'")
    phmm = MotifAwarePHMM(phmm_params)
    if require_anchored_repeat and bam_path is None and reads1 is not None:
        cache_records = select_informative_read_records(
            sample_id=sample_id,
            registry=registry,
            gaf_path=gaf_path,
            reads1=reads1,
            reads2=reads2,
            min_mapq=min_mapq,
        )
        locus_reads, informative_reads, _ = selected_reads_by_locus(cache_records)
        if use_read_cache and read_cache_dir is not None:
            write_read_cache(cache_records, read_cache_dir)
    else:
        locus_reads = collect_locus_reads(gaf_path, reads1, reads2, min_mapq=min_mapq, bam_path=bam_path)
        informative_reads = {}
    stats = {
        "likelihood_backend": likelihood_backend,
        "loci_with_reads": 0,
        "calls": 0,
        "mapped_reads": sum(len(v) for v in locus_reads.values()),
        "informative_reads": 0,
        "unmapped_or_unknown_locus_reads": 0,
        "novel_hypotheses": 0,
        "novel_accepted": 0,
    }
    novel_records: list[dict] = []
    with open_text(out_path, "wt") as out:
        writer = csv.DictWriter(out, fieldnames=OUTPUT_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for pointer_id in sorted(locus_reads):
            locus = registry.loci.get(pointer_id)
            if locus is None:
                continue
            mapped_reads = locus_reads[pointer_id]
            reads = (
                informative_reads.get(pointer_id, [])
                if require_anchored_repeat and informative_reads
                else (
                    [read for read in mapped_reads if has_anchored_repeat_evidence(read, locus)]
                    if require_anchored_repeat
                    else mapped_reads
                )
            )
            if not reads:
                continue
            stats["informative_reads"] += len(reads)
            stats["loci_with_reads"] += 1
            call, _, _ = call_locus(
                locus,
                reads,
                phmm,
                population_mix,
                theta,
                sigma,
                heterogeneity_threshold,
                prior_mode=prior_mode,
                use_length_smoothing=use_length_smoothing,
                likelihood_model=likelihood_backend,
                alignment_params=alignment_params,
                use_mixture_model=use_mixture_model,
            )
            final_alleles = locus.alleles
            if discover_novel:
                hypotheses, final_alleles = discover_novel_alleles(
                    locus,
                    reads,
                    phmm,
                    population_mix,
                    theta,
                    sigma,
                    call,
                    min_support=min_novel_support,
                    min_log_likelihood_gain=min_novel_log_gain,
                    min_gq=min_novel_gq,
                    likelihood_backend=likelihood_backend,
                    alignment_params=alignment_params,
                    prior_mode=prior_mode,
                    use_length_smoothing=use_length_smoothing,
                )
                for hypothesis in hypotheses:
                    novel_records.append({"pointer_id": pointer_id, **hypothesis.to_dict()})
                    stats["novel_hypotheses"] += 1
                    stats["novel_accepted"] += hypothesis.status == "accepted"
                if len(final_alleles) > len(locus.alleles):
                    call, _, _ = call_locus(
                        locus,
                        reads,
                        phmm,
                        population_mix,
                        theta,
                        sigma,
                        heterogeneity_threshold,
                        final_alleles,
                        prior_mode=prior_mode,
                        use_length_smoothing=use_length_smoothing,
                        likelihood_model=likelihood_backend,
                        alignment_params=alignment_params,
                        use_mixture_model=use_mixture_model,
                    )
            writer.writerow(
                {
                    "pointer_id": call.pointer_id,
                    "chrom": call.chrom,
                    "start": call.start,
                    "end": call.end,
                    "motif": call.motif,
                    "mapped_depth": len(mapped_reads),
                    "depth": call.depth,
                    "GT": call.gt,
                    "GQ": call.gq,
                    "allele1": call.allele1,
                    "allele2": call.allele2,
                    "posterior": f"{call.posterior:.10g}",
                    "log_likelihood": f"{call.log_likelihood:.10g}",
                    "log_prior": f"{call.log_prior:.10g}",
                    "log_score": f"{call.log_score:.10g}",
                    "log_mix": "NA" if call.log_mix is None else f"{call.log_mix:.10g}",
                    "D_PAPER": "NA" if call.d_paper is None else f"{call.d_paper:.10g}",
                    "D_LIKELIHOOD": "NA" if call.d_likelihood is None else f"{call.d_likelihood:.10g}",
                    "heterogeneity_threshold": "." if call.heterogeneity_threshold is None else call.heterogeneity_threshold,
                    "HET_FLAG": call.heterogeneity_flag,
                    "mixture_weights": (
                        "NA"
                        if call.mixture_weights is None
                        else json.dumps(call.mixture_weights, sort_keys=True)
                    ),
                }
            )
            stats["calls"] += 1
    if novel_out:
        with open_text(novel_out, "wt") as handle:
            for record in novel_records:
                handle.write(json.dumps(record, ensure_ascii=False, sort_keys=True) + "\n")
    return stats
