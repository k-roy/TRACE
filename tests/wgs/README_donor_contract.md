# Donor outcome contracts

Run `python -m pytest tests/wgs/test_donor_outcome.py tests/wgs/test_donor_outcome_contract.py` from the TRACE checkout with its dependencies installed. These tests author expected alleles and molecules independently and create their own temporary VCF/gVCF/FASTA inputs. They need no MAGESTIC or planning artifacts.

The samplesheet has seven required columns: `sample`, `chrom`, `donor_start_coord`, `donor_end_coord`, `wt_seq`, `design_donor`, `synthetic_donor`. Coordinates are 1-based inclusive on the actual reference. Optional `donor_strand` is `+` or `-` and applies together to all three supplied sequences.

A missing/blank strand permits WT/reference inference only when it uniquely determines orientation. If WT is its own reverse complement, an asymmetric design or empirical donor requires explicit strand. Fully symmetric WT/design/empirical sequences need no declaration because both orientations describe the same molecules. Contradictory or invalid declarations fail. The original samplesheet/specification retains the supplied strand; a returned genome-oriented copy uses `+` after flipping.

For example, reference WT `ACGT` and supplied design `CCGT` mean genomic `A1C` when strand is `+`, but genomic `T4G` when strand is `-`. Tests fix those expectations before providing VCFs. Strand is never inferred from observed variants. Existing upstream palindromic fixtures now explicitly declare the plus strand already specified by their original expected genomic alleles.
