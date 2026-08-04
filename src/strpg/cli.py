from __future__ import annotations

import argparse
import json
import sys
import warnings
from pathlib import Path

from .calibration import calibrate_threshold, load_threshold
from .defaults import DEFAULT_LIKELIHOOD_BACKEND
from .genotyper import genotype_dataset
from .graph import build_graph_and_registry
from .index import build_seed_index
from .mapper import map_bam, map_fastq
from .phmm import PHMMParameters
from .registry import AlleleRegistry
from .simulate import simulate_reads
from .update import update_registry


def _print_stats(stats: dict) -> None:
    print(json.dumps(stats, indent=2, sort_keys=True), file=sys.stderr)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="strpg",
        description="STR-PG pointer-graph short tandem repeat genotyping",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("build", help="Build fixed-topology pointer GFA and external allele registry")
    p.add_argument("--ref", required=True)
    p.add_argument("--str-catalog", required=True)
    p.add_argument("--out", required=True, help="Output GFA")
    p.add_argument("--registry", required=True, help="Output registry JSONL")
    p.add_argument("--vcf", default=None, help="Optional plain-text VCF for non-STR SNP/indel bubbles")
    p.add_argument("--build", default=None)
    p.add_argument("--flank", type=int, default=100)
    p.add_argument("--default-population", default="GLOBAL")
    p.add_argument("--initial-count-total", type=float, default=1000.0)

    p = sub.add_parser("index", help="Build flank-syncmer SQLite index")
    p.add_argument("--registry", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--k", type=int, default=15)
    p.add_argument("--s", type=int, default=5)
    p.add_argument("--t", type=int, default=2)

    p = sub.add_parser("map", help="Map FASTQ reads to pointer loci")
    p.add_argument("--index", required=True)
    source = p.add_mutually_exclusive_group(required=True)
    source.add_argument("--reads1")
    source.add_argument("--bam")
    p.add_argument("--reads2")
    p.add_argument("--out", required=True)
    p.add_argument("--delta", type=int, default=10)
    p.add_argument("--tau0", type=float, default=0.20)
    p.add_argument("--min-seeds", type=int, default=2)
    p.add_argument("--min-score-margin", type=float, default=0.25)
    p.add_argument("--max-seed-occurrences", type=int, default=200)
    p.add_argument("--threads", type=int, default=1, help="Compatibility option only; no parallel execution is implemented")
    p.add_argument(
        "--batch-size",
        "--batch_size",
        dest="batch_size",
        type=int,
        default=5000,
        help="Compatibility option only; batching is not implemented",
    )

    p = sub.add_parser(
        "genotype",
        help="Run Bayesian genotyping with SW (default) or experimental pHMM likelihoods",
        description=(
            "Genotype STR loci using a backend-independent informative-read set. "
            "sw: Default production backend based on Smith-Waterman alignment scores. "
            "phmm: Optional experimental motif-aware profile-HMM backend retained "
            "for reproducibility and comparative evaluation."
        ),
    )
    p.add_argument("--registry", required=True)
    p.add_argument("--gaf", required=True)
    source = p.add_mutually_exclusive_group(required=True)
    source.add_argument("--fq1")
    source.add_argument("--bam")
    p.add_argument("--fq2")
    p.add_argument("--out", required=True)
    p.add_argument("--popmix", default="GLOBAL=1")
    p.add_argument("--theta", type=float, default=0.01)
    p.add_argument("--sigma", type=float, default=10.0)
    p.add_argument(
        "--no-population-frequency",
        action="store_true",
        help="Use a uniform diploid genotype prior",
    )
    p.add_argument(
        "--no-balding-nichols",
        action="store_true",
        help="Use population-frequency HWE priors with theta fixed to zero",
    )
    p.add_argument(
        "--no-length-smoothing",
        action="store_true",
        help="Disable genotype-prior length smoothing",
    )
    p.add_argument(
        "--no-mixture-model",
        action="store_true",
        help="Disable optional mixture diagnostics without changing diploid GT/GQ",
    )
    p.add_argument("--min-mapq", type=int, default=0)
    p.add_argument(
        "--likelihood-backend",
        choices=("sw", "phmm"),
        default=DEFAULT_LIKELIHOOD_BACKEND,
        help="Allele likelihood backend (default: sw)",
    )
    p.add_argument(
        "--likelihood-model",
        choices=("sw", "phmm", "alignment_softmax"),
        default=None,
        help="Deprecated alias for --likelihood-backend",
    )
    p.add_argument("--mismatch-rate", type=float, default=0.01)
    p.add_argument("--flank-indel-open", type=float, default=0.003)
    p.add_argument("--repeat-indel-open", type=float, default=0.005)
    p.add_argument("--flank-indel-extend", type=float, default=0.10)
    p.add_argument("--repeat-indel-extend", type=float, default=0.30)
    p.add_argument(
        "--include-uninformative-flank-reads",
        action="store_true",
        help="Include pure-flank reads in allele likelihoods (not recommended; default uses flank-anchored repeat evidence only)",
    )
    p.add_argument(
        "--min-reference-log-gain",
        type=float,
        default=None,
        help="Deprecated and ignored; read selection no longer uses likelihood scores",
    )
    p.add_argument("--sample-id", default="sample")
    p.add_argument(
        "--read-cache-dir",
        help="Optional root for deterministic per-sample/per-locus read-cache TSV files",
    )
    p.add_argument(
        "--no-read-cache",
        action="store_true",
        help="Disable writing read cache files",
    )
    p.add_argument("--heterogeneity-threshold", type=float)
    p.add_argument("--heterogeneity-threshold-file")
    p.add_argument("--discover-novel", action="store_true")
    p.add_argument("--novel-out")
    p.add_argument("--min-novel-support", type=int, default=3)
    p.add_argument("--min-novel-log-gain", type=float, default=3.0)
    p.add_argument("--min-novel-gq", type=int, default=20)
    p.add_argument("--gfa", help="Accepted for compatibility; allele content is intentionally read from --registry, not GFA")
    p.add_argument("--freq", help="Legacy alias not used; migrate to --registry")
    p.add_argument("--region", help="Accepted for compatibility; filter the registry before running for production workflows")

    p = sub.add_parser("update", help="Update registry counts/frequencies and incorporate accepted novel alleles")
    p.add_argument("--registry-in", required=True)
    p.add_argument("--geno", required=True)
    p.add_argument("--registry-out", required=True)
    p.add_argument("--pop", required=True)
    p.add_argument("--min-GQ", type=int, default=20)
    p.add_argument("--base-total", type=float, default=1000.0)
    p.add_argument("--novel-hypotheses")
    p.add_argument("--pseudo-count-alpha", type=float, default=0.5)
    p.add_argument("--sample-id")
    p.add_argument("--allow-heterozygous", action="store_true")

    p = sub.add_parser("calibrate", help="Calibrate empirical heterogeneous-signal threshold from null calls")
    p.add_argument("--genotypes", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--column", default="D_PAPER")
    p.add_argument("--quantile", type=float, default=0.99)

    p = sub.add_parser("validate-registry", help="Check registry invariants")
    p.add_argument("--registry", required=True)

    p = sub.add_parser("simulate", help="Generate a small locus-specific FASTQ dataset")
    p.add_argument("--registry", required=True)
    p.add_argument("--pointer-id", required=True)
    p.add_argument("--weights", required=True, help="Example: 15=0.5,18=0.5 or 15=0.4,18=0.4,35=0.2")
    p.add_argument("--out1", required=True)
    p.add_argument("--out2")
    p.add_argument("--coverage", type=int, default=30)
    p.add_argument("--read-length", type=int, default=100)
    p.add_argument("--fragment-mean", type=int, default=180)
    p.add_argument("--fragment-sd", type=int, default=20)
    p.add_argument("--error-rate", type=float, default=0.002)
    p.add_argument("--stutter-rate", type=float, default=0.0)
    p.add_argument("--seed", type=int, default=13)
    return parser


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    if args.command == "build":
        build_graph_and_registry(args.ref, args.str_catalog, args.out, args.registry, args.build, args.flank, args.vcf, args.default_population, args.initial_count_total)
        return
    if args.command == "index":
        build_seed_index(args.registry, args.out, args.k, args.s, args.t)
        return
    if args.command == "map":
        if args.bam:
            _print_stats(map_bam(args.index, args.bam, args.out, args.delta, args.tau0, args.min_seeds, args.min_score_margin, args.max_seed_occurrences))
        else:
            _print_stats(map_fastq(args.index, args.reads1, args.out, args.reads2, args.delta, args.tau0, args.min_seeds, args.min_score_margin, args.max_seed_occurrences))
        return
    if args.command == "genotype":
        backend = args.likelihood_backend
        if args.likelihood_model is not None:
            alias = (
                "sw"
                if args.likelihood_model == "alignment_softmax"
                else args.likelihood_model
            )
            warnings.warn(
                "--likelihood-model is deprecated; use --likelihood-backend",
                FutureWarning,
                stacklevel=2,
            )
            if (
                args.likelihood_backend != DEFAULT_LIKELIHOOD_BACKEND
                and args.likelihood_backend != alias
            ):
                raise SystemExit(
                    "Conflicting --likelihood-backend and --likelihood-model values"
                )
            backend = alias
        if args.min_reference_log_gain is not None:
            warnings.warn(
                "--min-reference-log-gain is ignored because read selection is "
                "likelihood-backend-independent",
                FutureWarning,
                stacklevel=2,
            )
        threshold = args.heterogeneity_threshold
        if args.heterogeneity_threshold_file:
            threshold = load_threshold(args.heterogeneity_threshold_file)
        params = PHMMParameters(
            mismatch_rate=args.mismatch_rate,
            flank_indel_open=args.flank_indel_open,
            repeat_indel_open=args.repeat_indel_open,
            flank_indel_extend=args.flank_indel_extend,
            repeat_indel_extend=args.repeat_indel_extend,
        )
        prior_mode = (
            "uniform"
            if args.no_population_frequency
            else "hwe"
            if args.no_balding_nichols
            else "full"
        )
        effective_theta = 0.0 if prior_mode == "hwe" else args.theta
        _print_stats(
            genotype_dataset(
                registry_path=args.registry,
                gaf_path=args.gaf,
                reads1=args.fq1,
                out_path=args.out,
                reads2=args.fq2,
                population_mix_text=args.popmix,
                theta=effective_theta,
                sigma=args.sigma,
                min_mapq=args.min_mapq,
                phmm_params=params,
                heterogeneity_threshold=threshold,
                discover_novel=args.discover_novel,
                novel_out=args.novel_out,
                min_novel_support=args.min_novel_support,
                min_novel_log_gain=args.min_novel_log_gain,
                min_novel_gq=args.min_novel_gq,
                bam_path=args.bam,
                require_anchored_repeat=not args.include_uninformative_flank_reads,
                min_reference_log_gain=args.min_reference_log_gain,
                likelihood_backend=backend,
                sample_id=args.sample_id,
                read_cache_dir=args.read_cache_dir,
                use_read_cache=not args.no_read_cache,
                prior_mode=prior_mode,
                use_length_smoothing=not args.no_length_smoothing,
                use_mixture_model=not args.no_mixture_model,
            )        )
        return
    if args.command == "update":
        _print_stats(update_registry(args.registry_in, args.geno, args.registry_out, args.pop, args.min_GQ, args.base_total, args.novel_hypotheses, args.pseudo_count_alpha, args.sample_id, not args.allow_heterozygous))
        return
    if args.command == "calibrate":
        _print_stats(calibrate_threshold(args.genotypes, args.out, args.column, args.quantile))
        return
    if args.command == "validate-registry":
        errors = AlleleRegistry.load(args.registry).validate()
        if errors:
            print("\n".join(errors), file=sys.stderr)
            raise SystemExit(1)
        print("Registry validation passed.")
        return
    if args.command == "simulate":
        _print_stats(simulate_reads(args.registry, args.pointer_id, args.weights, args.out1, args.out2, args.coverage, args.read_length, args.fragment_mean, args.fragment_sd, args.error_rate, args.stutter_rate, args.seed))
        return
    raise AssertionError(args.command)


if __name__ == "__main__":
    main()
