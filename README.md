# STR-PG

STR-PG decouples STR allelic diversity from graph topology by representing
each STR locus with a fixed pointer node and storing allele sequences and
frequencies in an external registry.

This repository is the public reference implementation. Smith–Waterman (SW) is
the default production allele-likelihood backend. The corrected motif-aware
profile HMM (pHMM) is optional and experimental, retained for reproducibility
and comparative evaluation. Both consume the same backend-independent
informative-read set.

## Key design

- Fixed pointer nodes prevent graph topology from growing with the number of
  STR alleles.
- An external registry stores candidate sequences, repeat counts, provenance,
  and population frequencies.
- Syncmers and chaining assign reads to pointer loci.
- Read selection uses mapping, pairing, anchors, motif-compatible sequence, and
  fixed evidence rules—never an allele-likelihood score.
- A selected likelihood backend produces only a read-by-allele matrix.
- Diploid likelihoods, population priors, posterior/GQ, and optional mixture
  diagnostics are downstream and backend-agnostic.

## Likelihood-backend policy

`sw` is the default production backend based on bidirectional local-alignment
scores and a temperature-scaled log softmax. `phmm` is an optional experimental
motif-aware forward model. Current benchmark evidence did not support pHMM
superiority over SW.

```text
strpg genotype ... --likelihood-backend sw
strpg genotype ... --likelihood-backend phmm
```

Omitting the option is equivalent to `--likelihood-backend sw`.

## Installation

Python 3.9 or newer is required.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -U pip
python -m pip install -e ".[test]"
```

PowerShell:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -U pip
python -m pip install -e ".[test]"
```

BAM/CRAM input requires the optional `bio` extra.

## Quick start

```bash
strpg build \
  --ref examples/reference.fa \
  --str-catalog examples/str_catalog.tsv \
  --vcf examples/variants.vcf \
  --out demo_output/graph.gfa \
  --registry demo_output/registry.jsonl

strpg index \
  --registry demo_output/registry.jsonl \
  --out demo_output/index.sqlite

strpg map \
  --index demo_output/index.sqlite \
  --reads1 examples/demo_R1.fastq \
  --reads2 examples/demo_R2.fastq \
  --out demo_output/demo.gaf

strpg genotype \
  --registry demo_output/registry.jsonl \
  --gaf demo_output/demo.gaf \
  --fq1 examples/demo_R1.fastq \
  --fq2 examples/demo_R2.fastq \
  --sample-id demo \
  --read-cache-dir demo_output/read_cache \
  --out demo_output/genotypes.tsv
```

PowerShell users can run:

```powershell
PowerShell -ExecutionPolicy Bypass -File scripts\run_demo.ps1
PowerShell -ExecutionPolicy Bypass -File scripts\run_phmm_demo.ps1
```

## Full workflow

1. `build` creates a fixed-topology GFA and allele registry.
2. `index` builds the flank-syncmer SQLite index.
3. `map` assigns FASTQ or BAM/CRAM reads to pointer loci.
4. `genotype` selects informative reads independently of likelihood backend.
5. SW or pHMM computes the read-by-allele likelihood matrix.
6. The caller integrates diploid likelihoods and population-aware priors.
7. Optional mixture diagnostics describe heterogeneous allele signals without
   changing the diploid GT or GQ.
8. `update` can update registry counts from accepted calls.

## Input formats

- Reference: FASTA.
- STR catalog: tab-separated fields `build`, `chrom`, `start`, `end`, `motif`,
  `candidate_Ls`, and `pointer_id`.
- Variants: optional plain-text VCF for non-STR bubbles.
- Reads: paired or unpaired FASTQ; BAM/CRAM with the optional dependency.
- Mapping: GAF with pointer and mate tags produced by STR-PG.
- Registry: JSON Lines, one locus record per line.

See `docs/INPUT_FORMATS.md` and `docs/REGISTRY_FORMAT.md`.

## Output formats

The genotype TSV reports locus coordinates, mapped and informative depth, GT,
GQ, selected-genotype posterior, likelihood/prior/score terms, and optional
mixture diagnostics. Read-cache TSV files retain IDs, mates, sequence,
qualities, MAPQ, evidence flags, selection reasons, and source SHA256 values.
See `docs/OUTPUT_FORMATS.md`.

## SW backend

The production defaults are centralized in `strpg.defaults`:

- match: 2
- mismatch: -5
- gap: -5
- band: 500
- temperature: 12

For each read and allele, the maximum forward/reverse-complement local
alignment score is converted to a natural-log-compatible probability by a
temperature-scaled softmax. SW does not select reads, alter candidates, add a
prior, or enumerate genotypes.

## Optional pHMM backend

The pHMM uses forward M/I/D states, separate flank/repeat transition
parameters, and orientation marginalization. Its semi-global template-start
prior is applied once; termination does not divide by template length again.
It remains optional, experimental, and is not the default.

## Backend-independent read selection

Selection uses read-to-pointer mapping, MAPQ, paired-end relationships, flank
anchors, boundary-spanning motif evidence, sequence, and fixed evidence rules.
No SW score, pHMM likelihood, candidate ranking, posterior, prediction, or truth
is available to this stage. Both backends receive byte-identical cache content,
read order, and candidate order.

## Population-aware genotype prior

Registered population frequencies can be combined according to `--popmix`.
The genotype prior is applied once after all per-read diploid likelihoods are
summed. Length smoothing and Balding–Nichols-related controls are implemented
downstream of the likelihood matrix.

## Mixture diagnostics

Mixture fitting is an optional diagnostic. Disabling it changes only mixture
fields; it does not change diploid GT or GQ.

## Novel-allele handling

Novel candidates require conservative flank-supported repeat evidence.
Candidate creation is separate from the likelihood backend. Novel-allele
support is experimental and should be validated for the intended assay.

## Examples

The files under `examples/` are small deterministic synthetic fixtures and do
not contain private sequencing data. `scripts/run_demo.*` uses default SW;
`scripts/run_phmm_demo.*` explicitly requests pHMM.

## Tests

```bash
python -m compileall -q src
python -m pytest -q
python -m strpg.cli --help
python -m strpg.cli genotype --help
```

The release tests cover backend defaults and validation, read-selection
isolation, cache provenance and reproducibility, candidate ordering, pHMM
forward mathematics, genotype integration, private-path scanning, and file-size
limits.

## Reproducibility

Read caches are schema-versioned and deterministic. Source FASTQ and GAF
SHA256 values are stored with every cache record. Randomized synthetic
workflows use explicit seeds. See `docs/REPRODUCIBILITY.md`.

## Limitations

- Published validation currently concerns targeted loci and must not be
  interpreted as genome-wide validation.
- Pointer topology demonstrates a representation/storage design; it is not
  claimed by itself to improve genotype accuracy.
- pHMM is experimental and has not demonstrated superiority over SW.
- The optional BAM/CRAM path depends on `pysam`.
- A software license still requires copyright-holder approval; see
  `LICENSE_REQUIRED.md`.

## Citation

Use the author and title metadata in `CITATION.cff`. No DOI, journal acceptance,
or release DOI is asserted here.

## License

No license was present in the source project. Public distribution is blocked
until the copyright holders select one; see `LICENSE_REQUIRED.md`.

## Contributing

See `CONTRIBUTING.md`, `CODE_OF_CONDUCT.md`, and `SECURITY.md`.
