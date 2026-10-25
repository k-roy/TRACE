# Reconciliation functional tests

From the TRACE checkout, with its Python dependencies and pytest installed:

```sh
python -m pytest -q tests/wgs/test_reconcile_functional.py
```

The file creates and indexes its own small FASTA/BAM inputs with pysam. It has no
MAGESTIC, cluster, external-data, or audit-directory dependency. Expected alleles
and read molecules are specified before calling the implementation. Parameterized
cases cover exact SNP/MNV/indel identity, reference-equivalent indels, observed
repeat interruptions, CIGAR boundaries, quality and mapping filters, fragment
counting, incomplete designs, window contracts, and atomic CLI output.

The default locus policy requires MAPQ at least 20 and rejects MAPQ 255. Its
boundary tests are correctness checks of that explicit policy, not empirical
precision calibration. Legacy window maxima remain unassessed compatibility
evidence. Synthetic fixtures establish behavior on the declared molecules; they
do not establish accuracy on a biological cohort or completeness of a design
table. MAGESTIC owns tests of sample assignments and raw/donor-cleaned BAM wiring.
