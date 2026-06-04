#!/usr/bin/env python
"""
Parity harness: trace_crispr.wgs.sv_plasmid.clean_donor_reads (Half B) vs a
reference set of ``*.target_filtering.anno`` files produced by Shengdi Li's
original ``clean_donor_reads.pl``.

For each samplesheet row it streams the matching input BAM through
:func:`classify_donor_reads` (in ``--dup-mode parity`` to reproduce Shengdi's
fragment counts bit-for-bit) and compares the four ``.anno`` fields
``(variant_cov, frags_genome, frags_plasmid, frags_chimeric)`` against the
reference file. The ``.anno`` depends only on ``(input_bam, chrom, donor_start,
donor_end, blacklist)``, so results are cached per distinct key.

Optionally (``--check-bam-membership N``) it runs the full two-pass
:func:`clean_donor_reads` on the first ``N`` distinct samples and verifies that
the cleaned-BAM read-name membership matches the reference ``.clean.bam`` (this
is independent of ``--dup-mode``).

This is a manual integration test pointed at your own Sherlock WGS reference data
via CLI args (none of it ships with the package). Example::

    python parity_clean_donor_reads.py \\
        --samplesheet .../sample_BC_variant_fil.nns.v20250620.v2.txt \\
        --bam-dir     .../output3/bam \\
        --ref-anno-dir .../output3/clean_bam \\
        --blacklist   .../inputs/REDI_reads_mapped_regions.txt \\
        --samtools    /path/to/samtools \\
        [-n 0] [--check-bam-membership 10]
"""

import argparse
import os
import subprocess
import tempfile

from trace_crispr.wgs.sv_plasmid import (
    classify_donor_reads,
    clean_donor_reads,
    load_blacklist,
    parse_sam_line,
)

# samplesheet column names (== the v20250620.v2 header)
COL_SAMPLE = "sample"
COL_SHORT = "sample_name_short"
COL_CHROM = "chrom"
COL_DSTART = "donor_start_coord"
COL_DEND = "donor_end_coord"


def read_rows(samplesheet):
    with open(samplesheet) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {c: header.index(c) for c in (COL_SAMPLE, COL_SHORT, COL_CHROM, COL_DSTART, COL_DEND)}
        for line in fh:
            if not line.strip():
                continue
            c = line.rstrip("\n").split("\t")
            yield {
                "full": c[idx[COL_SAMPLE]],
                "short": c[idx[COL_SHORT]],
                "chrom": c[idx[COL_CHROM]],
                "dstart": int(c[idx[COL_DSTART]]),
                "dend": int(c[idx[COL_DEND]]),
            }


def read_ref_anno(path):
    """Return (variant_cov, frags_genome, frags_plasmid, frags_chimeric) from a .anno."""
    with open(path) as fh:
        fh.readline()  # header
        data = fh.readline().rstrip("\n").split("\t")
    return tuple(int(x) for x in data[:4])


def stream_records(bam_path, samtools):
    proc = subprocess.Popen([samtools, "view", "-h", bam_path], stdout=subprocess.PIPE, text=True)
    if proc.stdout is None:
        raise RuntimeError("failed to capture samtools stdout")
    try:
        for line in proc.stdout:
            if not line.startswith("@"):
                yield parse_sam_line(line)
    finally:
        proc.stdout.close()
        rc = proc.wait()
    if rc:
        raise subprocess.CalledProcessError(rc, proc.args)


def bam_read_names(bam_path, samtools):
    out = subprocess.run([samtools, "view", bam_path], stdout=subprocess.PIPE, text=True, check=True)
    return {ln.split("\t", 1)[0] for ln in out.stdout.splitlines()}


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--samplesheet", required=True)
    p.add_argument("--bam-dir", required=True)
    p.add_argument("--ref-anno-dir", required=True)
    p.add_argument("--blacklist", required=True)
    p.add_argument("--samtools", default="samtools")
    p.add_argument("--dup-mode", choices=("parity", "standard"), default="parity")
    p.add_argument("-n", type=int, default=0, help="max testable rows to check (0 = all)")
    p.add_argument("--check-bam-membership", type=int, default=0, metavar="N",
                   help="also verify cleaned-BAM membership on the first N distinct samples")
    p.add_argument("--show", type=int, default=20, help="max mismatches to print")
    args = p.parse_args(argv)

    blacklist = load_blacklist(args.blacklist)

    cache = {}  # (short, chrom, dstart, dend) -> computed (vcov,G,P,C)
    n_tested = 0
    n_exact = 0
    mismatches = []
    membership_done = 0
    membership_ok = 0
    membership_fail = []

    for row in read_rows(args.samplesheet):
        bam = os.path.join(args.bam_dir, row["short"] + ".sort.dedup.recal2.bam")
        ref_anno = os.path.join(args.ref_anno_dir, row["full"] + ".clean.bam.target_filtering.anno")
        if not os.path.exists(bam) or not os.path.exists(ref_anno):
            continue

        key = (row["short"], row["chrom"], row["dstart"], row["dend"])
        if key in cache:
            mine = cache[key]
        else:
            res = classify_donor_reads(
                stream_records(bam, args.samtools),
                row["chrom"], row["dstart"], row["dend"], blacklist,
                dup_mode=args.dup_mode,
            )
            mine = (res.variant_cov, res.frags_genome, res.frags_plasmid, res.frags_chimeric)
            cache[key] = mine

            if args.check_bam_membership and membership_done < args.check_bam_membership:
                with tempfile.TemporaryDirectory() as td:
                    out_bam = os.path.join(td, "clean.bam")
                    clean_donor_reads(bam, out_bam, row["chrom"], row["dstart"], row["dend"],
                                      args.blacklist, dup_mode=args.dup_mode, samtools=args.samtools)
                    ref_clean = os.path.join(args.ref_anno_dir, row["full"] + ".clean.bam")
                    if os.path.exists(ref_clean):
                        mine_names = bam_read_names(out_bam, args.samtools)
                        ref_names = bam_read_names(ref_clean, args.samtools)
                        membership_done += 1
                        if mine_names == ref_names:
                            membership_ok += 1
                        else:
                            membership_fail.append(
                                (row["full"], len(mine_names ^ ref_names),
                                 sorted(mine_names - ref_names)[:3], sorted(ref_names - mine_names)[:3])
                            )

        ref = read_ref_anno(ref_anno)
        n_tested += 1
        if mine == ref:
            n_exact += 1
        else:
            mismatches.append((row["full"], mine, ref))

        if args.n and n_tested >= args.n:
            break

    print(f"distinct (bam,coords) computed : {len(cache)}")
    print(f"anno comparisons               : {n_tested}")
    print(f"exact matches                  : {n_exact}")
    print(f"mismatches                     : {len(mismatches)}")
    if mismatches:
        print("\n-- mismatches (full_sample | mine(vcov,G,P,C) | ref) --")
        for full, mine, ref in mismatches[: args.show]:
            print(f"  {full}\n      mine={mine}  ref={ref}")
    if args.check_bam_membership:
        print(f"\nBAM membership checked         : {membership_done}")
        print(f"BAM membership OK              : {membership_ok}")
        if membership_fail:
            print("-- membership failures (full | n_diff | mine_only[:3] | ref_only[:3]) --")
            for full, ndiff, mo, ro in membership_fail[: args.show]:
                print(f"  {full} | {ndiff} | {mo} | {ro}")

    ok = (len(mismatches) == 0) and (not args.check_bam_membership or not membership_fail)
    print("\nPARITY:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
