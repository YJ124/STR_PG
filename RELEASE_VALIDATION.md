# Release validation procedure

Run from the repository root in a clean Python 3.9+ environment:

```text
python -m pip install -e ".[test]"
python -m compileall -q src
python -m pytest -q
python -m strpg.cli --help
python -m strpg.cli genotype --help
```

Then run the default SW and explicit pHMM synthetic demos:

```text
PowerShell -ExecutionPolicy Bypass -File scripts\run_demo.ps1
PowerShell -ExecutionPolicy Bypass -File scripts\run_phmm_demo.ps1
```

Required checks:

1. both commands exit zero;
2. SW statistics record `"likelihood_backend": "sw"`;
3. pHMM statistics record `"likelihood_backend": "phmm"`;
4. both caches have the same SHA256;
5. both genotype TSV files have one valid data row;
6. no pHMM performance advantage is required;
7. repository scans find no private path, credential, partial file, cache, or
   large non-public data;
8. documentation links resolve;
9. `pyproject.toml`, `strpg.__version__`, `CITATION.cff`, and CHANGELOG agree;
10. copyright holders approve a LICENSE before public upload.

GitHub Actions repeats compilation, tests, CLI help, and both demos on Ubuntu
and Windows with Python 3.9 and 3.12. Local success does not claim hosted CI
success until the workflow has actually run on GitHub.
