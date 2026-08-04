# Architecture

STR-PG represents each STR locus with one stable pointer node. Allele sequences
and frequencies live in an external JSONL registry, so adding an allele does
not require adding an STR path to graph topology.

The runtime stages are:

1. reference/catalog parsing and pointer-GFA construction;
2. registry construction and validation;
3. flank-syncmer indexing;
4. seed lookup and chaining to assign reads to pointer loci;
5. backend-independent read selection and immutable cache creation;
6. SW or experimental pHMM read-by-allele likelihood calculation;
7. diploid genotype enumeration and one-time prior integration;
8. posterior/GQ calculation;
9. optional mixture diagnostics.

Candidate order is inherited from the registry and preserved through the
likelihood matrix and genotype enumeration. Likelihood backends cannot select
reads, change candidates, add priors, or inspect truth. Runtime code under
`src/strpg` has no dependency on experimental directories.
