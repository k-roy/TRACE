# bc1 boundary regression tests

`test_bc1_boundaries.py` checks molecular barcode recovery independently of
alignment through a degenerate reference N20. Run from the TRACE checkout:

```sh
python -m pytest tests/wgs/test_bc1_boundaries.py -q
```

The two configured profiles are `GGTCATGCTC/N20/ATAACAGGGT` (pF746-style)
and `GGTCATGCTC/N20/CCCTAGGATA` (pF835-style), with their reverse complements.
These are examples of per-track configuration, not universal era assignments.

Four frozen real reads in `fixtures/bc1_real_reads.json` contain a complete exact
left10/barcode20/right10 cassette. Each expected barcode was fixed from that
sequence before re-executing either caller. The tests reconstruct indexed BAMs
with the original sequence, CIGAR, flags, qualities and pair coordinates, then
exercise `call_sample`. Source records came from the MAGESTIC v11 audit; local
Oak paths are provenance only and are not needed to run these fixtures.

Other tests cover both profiles and orientations, both strands, recoverable
soft clips next to hard clips, unavailable hard-clipped bases, partial and
incorrect flanks, conflicting widths, ambiguous frames and a sequence motif
without a locus anchor. A deterministic 500-example test generates molecular
barcodes first and independently varies only the CIGAR through the variable
region. These establish decoding contracts, not well-level biological accuracy.
