# Contributing

Contributions are welcome through GitHub issues and pull requests after the
repository is published.

1. Create an isolated environment.
2. Install with `python -m pip install -e ".[test]"`.
3. Add tests for every behavior change.
4. Run `python -m compileall -q src` and `python -m pytest -q`.
5. Keep runtime code independent of local paths, truth genotypes, and
   `experiments/`.

Changes to likelihood mathematics, read selection, candidate ordering, or
priors must include focused unit tests and a reproducibility note. Do not add
private sequencing data or credentials.
