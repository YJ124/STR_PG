# Backend-independent read selection

Selection is implemented in `strpg.read_selection`. Inputs are mapping
assignments, MAPQ, pairing, read sequence/quality, locus flanks, and motif.
A read is informative when fixed flank/repeat boundary evidence is observed.
When one mate carries boundary evidence, both consistently mapped mates are
retained by the paired interval-union rule.

Selection does not import or instantiate a likelihood backend and cannot use
allele scores, candidate ranking, priors, posterior, predicted genotype, or
truth. Records are sorted by locus, read ID, and mate before cache writing.

The version 1.0 TSV cache records schema version, sample/locus/read IDs, mate,
sequence, base quality, MAPQ, selection status/reason, anchor/repeat flags,
sequence SHA256, source FASTQ SHA256, and source GAF SHA256. Cache writing can
be disabled with `--no-read-cache`.
