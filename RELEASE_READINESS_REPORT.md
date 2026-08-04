# STR-PG 2.0.0 release-readiness report

Date: 2026-07-29

## Outcome

The independent `STR-PG-sw` repository candidate is code-complete and locally
validated. It is technically suitable for GitHub after one mandatory human
decision: the copyright holders must approve a software license. Hosted GitHub
Actions has not run because this work intentionally did not create or push a
repository.

## Repository construction

The candidate was created as a new directory without modifying the source
project. Copied content was restricted to production Python modules, existing
core tests, pHMM mathematical tests, small synthetic examples, CLI compatibility
wrappers, requirements, citation metadata, and repository attributes.

Not copied: experimental result trees, archives, raw FASTQ/GAF collections,
SQLite indexes, read caches, partial files, local environment caches, IDE
configuration, credentials, or private data. Demo outputs and interpreter
caches created during validation were removed from the candidate afterward.

## Implemented architecture

- Unified interface:
  `src/strpg/likelihoods/base.py`
- Default production SW backend:
  `src/strpg/likelihoods/sw.py`
- Optional experimental pHMM backend:
  `src/strpg/likelihoods/phmm.py`
- Backend factory:
  `src/strpg/likelihoods/__init__.py`
- Central defaults:
  `src/strpg/defaults.py` and `config/defaults.yaml`
- Backend-independent selection/cache:
  `src/strpg/read_selection.py`
- Genotype integration:
  `src/strpg/genotyper.py`
- CLI policy and compatibility:
  `src/strpg/cli.py`

SW and pHMM receive the same ordered selected sequences and the same ordered
candidate list. Backends only create a natural-log-compatible read-by-allele
matrix. Genotype enumeration, one-time prior integration, posterior/GQ, and
mixture diagnostics remain downstream.

## Defaults

The formal defaults are:

- likelihood backend: SW
- backend-independent read selection: enabled
- population frequency: enabled
- Balding–Nichols correction: enabled
- length smoothing: enabled
- mixture diagnostics: enabled
- SW: match 2, mismatch -5, gap -5, band 500, temperature 12
- pHMM: audited defaults 0.01 mismatch, 0.003/0.10 flank indel
  open/extension, 0.005/0.30 repeat indel open/extension

The non-generalizing calibration P05 was not made a production default.

## Read-cache schema

Schema version 1.0 records sample ID, locus ID, read ID, mate, sequence, base
quality, MAPQ, selection flag/reason, left/right anchors, repeat evidence,
sequence SHA256, source FASTQ SHA256, and source GAF SHA256. Output is
deterministically sorted and byte-identical for identical inputs. Cache writing
can be disabled.

## CLI and compatibility

Current production use:

```text
strpg genotype ... --likelihood-backend sw
```

Experimental pHMM:

```text
strpg genotype ... --likelihood-backend phmm
```

Omitting the option selects SW. The old `--likelihood-model` option remains a
deprecated one-cycle alias and emits a warning. `alignment_softmax` maps to
`sw`; conflicts with an explicit new option are rejected.

Version is 2.0.0 because changing the default backend changes inference
behavior. `pyproject.toml`, package version, citation metadata, and CHANGELOG
agree.

## Local validation

- Editable installation: passed
- `compileall`: passed
- Pytest: **38 passed**
- CLI top-level help: passed
- Genotype help: passed and reports default SW
- Default SW synthetic demo: exit 0, one call, GT 10/12,
  backend recorded as `sw`
- Explicit pHMM synthetic demo: exit 0, one call, GT 10/12,
  backend recorded as `phmm`
- SW/pHMM read-cache SHA256 equality: passed
- SW/pHMM candidate order and matrix shape: passed
- Corrected pHMM brute-force/boundary tests: passed
- Mixture-disabled GT/GQ invariance: passed

The pHMM demo is not claimed to outperform SW.

## CI status

`.github/workflows/tests.yml` defines Ubuntu and Windows jobs for Python 3.9
and 3.12. It installs the test extra, compiles source, runs pytest, checks CLI
help, and runs both demos without large downloads.

Status: **configured but not hosted-run**. A green GitHub badge or hosted CI
success must not be claimed until the workflow runs after upload.

## Documentation

- `README.md`
- `README.zh-CN.md`
- `docs/ARCHITECTURE.md`
- `docs/LIKELIHOOD_BACKENDS.md`
- `docs/READ_SELECTION.md`
- `docs/REGISTRY_FORMAT.md`
- `docs/INPUT_FORMATS.md`
- `docs/OUTPUT_FORMATS.md`
- `docs/REPRODUCIBILITY.md`
- `docs/MIGRATION_TO_SW_DEFAULT.md`
- `docs/PHMM_AUDIT_SUMMARY.md`
- `RELEASE_VALIDATION.md`

Documentation uses relative examples and makes no pHMM superiority,
pointer-accuracy, genome-wide validation, journal acceptance, or DOI claim.

## Repository hygiene

- Private absolute paths found: 0
- Credential/key signatures found: 0
- Files larger than 2 MiB found: 0
- GAF/BAM/CRAM/SQLite/partial release files found: 0
- Python/pytest caches retained: 0
- Batch experiment results or archives retained: 0
- Small synthetic FASTQ fixtures: intentionally retained

## License

The source project had no LICENSE. No legal license was invented. The mandatory
release blocker is documented in `LICENSE_REQUIRED.md`.

## Scientific blockers

The corrected pHMM remains experimental, and existing evidence does not support
superiority over SW. This does not block the SW-default software architecture,
but it blocks any pHMM-advantage claim. Existing targeted-locus evidence must
not be described as genome-wide validation. Broader independent biological
validation remains future work.

## Release decision

**Technical GitHub readiness: YES.**

**Unconditional public-release readiness: NO, pending license approval and the
first hosted CI run.**

Before upload, the user must:

1. approve and add a software LICENSE;
2. review author/citation metadata;
3. inspect the release manifests;
4. initialize the repository and run hosted CI;
5. confirm that no additional unpublished material should be included.

No `git push`, repository creation, GitHub release, or PyPI upload was executed.
