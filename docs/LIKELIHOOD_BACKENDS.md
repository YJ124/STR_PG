# Likelihood backends

Both backends implement `AlleleLikelihoodBackend.score_reads` and return a
read-by-allele matrix of natural-log-compatible values where larger is better.

## SW

SW is the default production backend. It computes local-alignment scores for
the read and reverse complement, retains the larger score, and applies a
temperature-scaled log softmax across candidates. Defaults are match 2,
mismatch -5, gap -5, band 500, and temperature 12.

## pHMM

pHMM is optional and experimental. It uses forward M/I/D states with
flank/repeat-specific transitions and orientation marginalization. The
template-start prior is applied exactly once; termination sums mutually
exclusive end states without another template-length division.

pHMM is retained for reproducibility and comparative evaluation. Current
benchmark evidence did not support superiority over SW.
