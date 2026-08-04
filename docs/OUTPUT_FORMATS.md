# Output formats

The genotype TSV contains pointer ID, coordinates, motif, mapped depth,
informative depth, GT, GQ, alleles, selected-genotype posterior, log
likelihood, log prior, log score, mixture likelihood, heterogeneity statistics,
and mixture weights.

Likelihood and prior terms are natural-log compatible. The prior is added once
after per-read diploid likelihoods are summed.

Read cache files are written to:

```text
read_cache/<sample_id>/<locus_id>.tsv
```

Their schema is documented in `READ_SELECTION.md`. Novel hypotheses, when
requested, are JSON Lines and remain provisional unless all acceptance
criteria pass.
