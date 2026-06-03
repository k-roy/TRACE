#!/usr/bin/env python
"""
Parity harness: trace_crispr.wgs.donor_outcome vs a reference precalling output.

Runs the Python pre-caller over a subset (or all) of the samples in a samplesheet
and compares, field-for-field, against a reference `on_target_precalling.txt`
(produced by the original Perl ``precalling_target_outcome.pl``). Reports any
mismatches with full detail. Intended as a manual integration test, pointed at
your own WGS reference data via CLI args (it needs the per-sample VCF/gVCF + the
genome FASTA + the reference summary, none of which ship with the package).

Usage:
    python parity_precalling.py \\
        --samplesheet sample_BC_variant.tsv \\
        --workdir results_dir \\            # contains vcf/{s}.fil.vcf + gvcf/{s}.g.vcf
        --genome combined_reference.fa \\
        --refout on_target_precalling.txt \\
        [-n 150]                            # samples to check (0 = all)
"""

import argparse
import sys
from pathlib import Path

from trace_crispr.wgs.donor_outcome import (
    load_fasta,
    precall_sample,
    read_samplesheet,
)


def load_reference(path):
    ref = {}
    with open(path) as fh:
        fh.readline()  # header
        for line in fh:
            f = line.rstrip("\n").split("\t")
            # pad to 5 fields (trailing empty cols may be dropped by split on rstrip)
            while len(f) < 5:
                f.append("")
            ref[f[0]] = (f[1], f[2], f[3], f[4])
    return ref


def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--samplesheet", required=True, type=Path,
                   help="TSV with sample/chrom/donor_start_coord/donor_end_coord/"
                        "wt_seq/design_donor/synthetic_donor columns")
    p.add_argument("--workdir", required=True, type=Path,
                   help="dir containing vcf/{sample}.fil.vcf and gvcf/{sample}.g.vcf")
    p.add_argument("--genome", required=True, type=Path, help="reference genome FASTA")
    p.add_argument("--refout", required=True, type=Path,
                   help="reference on_target_precalling.txt to compare against")
    p.add_argument("-n", "--num-samples", type=int, default=150,
                   help="number of samples to check (0 = all; default 150)")
    args = p.parse_args()

    n = args.num_samples
    ref = load_reference(args.refout)
    refs_fasta = load_fasta(str(args.genome))
    specs = read_samplesheet(args.samplesheet)

    checked = mism = missing = 0
    mismatches = []
    for spec in specs:
        if spec.sample not in ref:
            continue
        vcf = args.workdir / "vcf" / f"{spec.sample}.fil.vcf"
        gvcf = args.workdir / "gvcf" / f"{spec.sample}.g.vcf"
        if not vcf.exists() or not gvcf.exists():
            missing += 1
            continue
        try:
            out = precall_sample(spec, vcf, gvcf, refs_fasta[spec.chrom])
        except Exception as e:  # noqa: BLE001
            mismatches.append((spec.sample, "EXCEPTION", str(e)))
            mism += 1
            checked += 1
            continue
        got = (
            out.editing_outcome,
            ",".join(out.designed_variants),
            ",".join(out.synthetic_errors),
            ",".join(out.spontaneous_variants),
        )
        exp = ref[spec.sample]
        checked += 1
        if got != exp:
            mism += 1
            if len(mismatches) < 25:
                cols = ["outcome", "designed", "synthetic", "spontaneous"]
                diffs = [
                    f"    {cols[i]}:\n      exp={exp[i]!r}\n      got={got[i]!r}"
                    for i in range(4)
                    if got[i] != exp[i]
                ]
                mismatches.append((spec.sample, "MISMATCH", "\n".join(diffs)))
        if n and checked >= n:
            break

    print(f"checked={checked}  matched={checked - mism}  mismatched={mism}  "
          f"missing_inputs={missing}")
    for sample, kind, detail in mismatches:
        print(f"\n[{kind}] {sample}\n{detail}")
    print(
        f"\nPARITY: {'PASS' if mism == 0 and checked > 0 else 'FAIL'} "
        f"({checked - mism}/{checked})"
    )
    return 1 if mism else 0


if __name__ == "__main__":
    raise SystemExit(main())
