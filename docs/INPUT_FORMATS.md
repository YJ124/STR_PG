# Input formats

- **FASTA:** standard named reference sequences.
- **STR catalog TSV:** `build`, `chrom`, `start`, `end`, `motif`,
  `candidate_Ls`, `pointer_id`; coordinates are one-based inclusive.
- **VCF:** optional plain-text non-STR SNP/indel records for graph bubbles.
- **FASTQ:** four-line records; paired files must contain matching IDs.
- **BAM/CRAM:** optional and requires `pysam`; secondary/supplementary records
  are ignored.
- **GAF:** at least 12 fields; STR-PG uses `pi`, `as`, and `mt` tags when
  present.
- **Registry:** JSONL described in `REGISTRY_FORMAT.md`.

Compressed text inputs are supported where `strpg.utils.open_text` recognizes
the suffix.
