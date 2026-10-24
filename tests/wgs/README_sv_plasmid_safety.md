# SV and plasmid safety regressions

Run from the TRACE repository with its WGS dependencies installed:

```sh
python -m pytest tests/wgs/test_sv_plasmid_safety.py -q
```

These 83 parametrized cases require `pytest`, `pysam`, and `samtools` on `PATH`. All FASTA, BAM, BAI/CSI and blacklist fixtures are generated in pytest temporary directories. No private data, audit workspace, download or MAGESTIC installation is required. Two bounded threads exercise concurrent writers; other tests are serial.

Expected sequences and read geometry are specified in the fixtures before the caller executes. Coverage includes two-sided junction anchors; CIGAR `D`/`N`/`=`/`X`; MAPQ, primary status and mate validity; exact blacklist/reference identity; standard versus explicit legacy parity counting; full cleaned-BAM preservation and indexed random access; deterministic output bytes; failed conversion/indexing; same-output writer serialization; and reference construction. The sparse CSI case uses a long header and one short read, not a large allocated genome.

This establishes synthetic functional contracts. It does not measure biological integration sensitivity, empirical MAPQ calibration, or advisory-lock behavior on network filesystems. Writers must honor the same POSIX lock protocol; the persistent empty `.lock` file is not a completion marker.
