# Reproducibility

- Install from a recorded commit and record `strpg.__version__`.
- Preserve registry, FASTQ/BAM, and GAF SHA256 values.
- Use the versioned read cache when comparing likelihood backends.
- Keep candidate ordering unchanged.
- Record CLI arguments, population mixture, theta, sigma, and backend.
- Use explicit random seeds for simulation.
- Do not select parameters using held-out evaluation samples.

Repeated cache construction from identical inputs is byte-identical. The public
demo is deterministic and contains only synthetic data. CI runs offline against
these small fixtures.
