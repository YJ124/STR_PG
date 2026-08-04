# Registry format

The registry is JSON Lines with one locus object per line. Required locus
content includes pointer ID, reference build, chromosome, one-based inclusive
STR coordinates, motif, left/right flanks, anchor length, and ordered alleles.

Each allele contains a stable ID, repeat count, complete sequence, repeat
sequence, motif, zero-based half-open repeat coordinates within the allele
template, population frequencies/counts, status, source, and optional metadata.

The complete allele sequence invariant is:

```text
left_flank + motif repeated repeat_count times + right_flank
```

Candidate order is semantically significant and is preserved during scoring.
Use `strpg validate-registry --registry PATH` to check invariants.
