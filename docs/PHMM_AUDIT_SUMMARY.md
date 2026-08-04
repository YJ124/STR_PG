# pHMM audit summary

The former semi-global forward implementation applied a uniform
template-start probability and then divided the summed terminal probability by
template length a second time. Brute-force path enumeration demonstrated the
extra template-length factor.

The minimal correction retains start uniformization and removes the second
termination division. Brute-force, probability-normalization, transition,
boundary, orientation, candidate-index, and diploid-integration tests cover
the corrected implementation.

pHMM remains optional and experimental for reproducibility. Independent audit
benchmarks did not support a claim that it outperforms SW, so SW is the
production default.
