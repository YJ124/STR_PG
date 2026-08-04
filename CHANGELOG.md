# Changelog

## 2.0.0 — 2026-07-29

- Changed the default genotype likelihood backend from pHMM to
  Smith–Waterman (`sw`). This affects default inference behavior and therefore
  warrants a major version.
- Added a unified allele-likelihood backend interface.
- Retained the corrected motif-aware pHMM as an explicit experimental backend.
- Removed likelihood-dependent informative-read selection.
- Added deterministic, versioned read-cache TSV files with read IDs, mates,
  qualities, mapping quality, evidence reasons, and source hashes.
- Added CLI conflict detection and a deprecated `--likelihood-model` alias.
- Added public synthetic examples, release tests, documentation, and CI.

## 1.0.0

- Initial manuscript-aligned reference implementation.
