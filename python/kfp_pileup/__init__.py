"""kfp-pileup: read recruitment and crude pileups downstream of a kfp profile (plan, phase 13).

A separate tool from kfp (plan: Tools and scope). It reads kfp's profile, index tables and
target sequences (``kmer_functional_profiler.store``), builds a per-sample reference
(:mod:`kfp_pileup.reference`), leaves alignment to DIAMOND (amino-acid indexes) or bwa-mem2
(nucleotide indexes), and turns the alignments into read posteriors, pileups and a per-unit
identity table (:mod:`kfp_pileup.pileup`). Crude by design: enough to evaluate kfp as a
front end, not a variant caller.
"""
