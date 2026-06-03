#!/usr/bin/env python
"""
Parity harness: trace_crispr.wgs.donor_outcome vs Shengdi's reference output.

Runs the Python pre-caller over a subset (or all) of the samples in Shengdi's
samplesheet and compares, field-for-field, against his reference
`on_target_precalling.txt`. Reports any mismatches with full detail.

Usage:
    python parity_check.py [N]        # N samples (default 150; 0 = all)
"""

import sys
from pathlib import Path

from trace_crispr.wgs.donor_outcome import (
    load_fasta,
    precall_sample,
    read_samplesheet,
)

SH = Path("/path/to/projects")
OM = SH / "WGS_MAGESTIC_QTL_outcome_mapping"
BC = SH / "WGS_MAGESTIC_QTL_bc_calling"
SAMPLESHEET = OM / "inputs" / "sample_BC_variant_fil.nns.v20250620.v2.txt"
WORKDIR = BC / "output3"
GENOME = OM / "inputs" / "fasta" / "combined_reference.v2.fa"
REFOUT = BC / "output3" / "summaries" / "on_target_precalling.txt"


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
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 150
    ref = load_reference(REFOUT)
    refs_fasta = load_fasta(str(GENOME))
    specs = read_samplesheet(SAMPLESHEET)

    checked = mism = missing = 0
    mismatches = []
    for spec in specs:
        if spec.sample not in ref:
            continue
        vcf = WORKDIR / "vcf" / f"{spec.sample}.fil.vcf"
        gvcf = WORKDIR / "gvcf" / f"{spec.sample}.g.vcf"
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

    print(
        f"checked={checked}  matched={checked - mism}  mismatched={mism}  missing_inputs={missing}"
    )
    for sample, kind, detail in mismatches:
        print(f"\n[{kind}] {sample}\n{detail}")
    print(
        f"\nPARITY: {'PASS' if mism == 0 and checked > 0 else 'FAIL'} ({checked - mism}/{checked})"
    )
    return 1 if mism else 0


if __name__ == "__main__":
    raise SystemExit(main())
