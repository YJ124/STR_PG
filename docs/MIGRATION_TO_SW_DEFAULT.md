# Migration to the SW default

Version 1.x used pHMM as the implicit genotype likelihood. Version 2.0 changes
the default to SW, so this is a behavior-changing major release.

- Current default: `strpg genotype ...`
- Explicit current default: `strpg genotype ... --likelihood-backend sw`
- Reproduce the corrected old backend:
  `strpg genotype ... --likelihood-backend phmm`

The deprecated `--likelihood-model` alias remains for one compatibility cycle,
emits a warning, and maps `alignment_softmax` to `sw`. Conflicting new and old
options are rejected.

Genotype output columns are unchanged; run statistics now record
`likelihood_backend`. Historical pHMM experiments do not represent the current
default model.
