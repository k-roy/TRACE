#!/usr/bin/env python3
"""
Align unique R1/R2 pairs to reference using BWA-MEM paired-end mode.

NEW APPROACH: Create separate R1.fastq and R2.fastq files and run
BWA-MEM in paired-end mode to properly handle:
- Insertions at cut site (donor capture)
- Reads that can't merge
- Paired-end alignment benefits (insert size, proper pairing flags)

Supports three modes:
- Paired-end: Both R1 and R2 present
- Single-end: Only R1 present (R2 is empty string)

The seq_hash is stored in the read name for later lookup.
"""

import gzip
import subprocess
import tempfile
from pathlib import Path

import pandas as pd


def main():
    unique_seqs_path = snakemake.input.unique_seqs
    ref_path = snakemake.input.ref
    output_bam = snakemake.output.bam
    output_bai = snakemake.output.bai
    log_path = Path(str(snakemake.log))
    threads = snakemake.threads

    log_path.parent.mkdir(parents=True, exist_ok=True)
    Path(output_bam).parent.mkdir(parents=True, exist_ok=True)

    # Load unique R1/R2 pairs
    df = pd.read_csv(unique_seqs_path, sep='\t', compression='gzip')

    # Detect sequencing mode
    # For paired-end mode, ALL reads must have R2 (BWA requires matched pairs)
    # If mixed (some have R2, some don't), fall back to single-end mode
    def has_valid_r2(val):
        return val and not pd.isna(val) and str(val).strip() != ''

    r2_status = [has_valid_r2(row['r2_sequence']) for _, row in df.iterrows()]
    all_have_r2 = all(r2_status)
    any_have_r2 = any(r2_status)

    if all_have_r2:
        mode = "paired-end"
        use_r2 = True
    elif any_have_r2:
        mode = "mixed (using single-end)"
        use_r2 = False  # Can't use mixed with BWA paired-end
    else:
        mode = "single-end"
        use_r2 = False

    # Create temporary FASTQ files
    # Read name = seq_hash for easy lookup (same for R1 and R2)
    r1_fq = tempfile.NamedTemporaryFile(mode='w', suffix='_R1.fq', delete=False)
    r1_fq_path = r1_fq.name

    if use_r2:
        r2_fq = tempfile.NamedTemporaryFile(mode='w', suffix='_R2.fq', delete=False)
        r2_fq_path = r2_fq.name
    else:
        r2_fq_path = None

    # Write FASTQ records
    for _, row in df.iterrows():
        r1_seq = row['r1_sequence']
        seq_hash = row['seq_hash']

        # Write R1 record
        r1_fq.write(f"@{seq_hash}\n")
        r1_fq.write(f"{r1_seq}\n")
        r1_fq.write("+\n")
        r1_fq.write("I" * len(r1_seq) + "\n")  # Placeholder quality

        # Write R2 record only in pure paired-end mode
        if use_r2:
            r2_seq = row['r2_sequence']
            r2_fq.write(f"@{seq_hash}\n")
            r2_fq.write(f"{r2_seq}\n")
            r2_fq.write("+\n")
            r2_fq.write("I" * len(r2_seq) + "\n")  # Placeholder quality

    r1_fq.close()
    if use_r2:
        r2_fq.close()

    # Run BWA alignment in appropriate mode
    if use_r2:
        # Paired-end mode
        cmd = f"bwa mem -t {threads} {ref_path} {r1_fq_path} {r2_fq_path} | samtools sort -@ 2 -o {output_bam} -"
    else:
        # Single-end mode
        cmd = f"bwa mem -t {threads} {ref_path} {r1_fq_path} | samtools sort -@ 2 -o {output_bam} -"

    with open(log_path, 'w') as log:
        log.write(f"Aligning {len(df)} unique sequences\n")
        log.write(f"Mode: {mode}\n")
        log.write(f"R2 stats: {sum(r2_status)}/{len(r2_status)} reads have R2\n")
        log.write(f"Command: {cmd}\n\n")

        result = subprocess.run(
            cmd,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )

        log.write(f"STDERR:\n{result.stderr.decode()}\n")

        if result.returncode != 0:
            raise RuntimeError(f"BWA alignment failed: {result.stderr.decode()}")

    # Index BAM
    subprocess.run(f"samtools index {output_bam}", shell=True, check=True)

    # Clean up temp files
    Path(r1_fq_path).unlink()
    if r2_fq_path:
        Path(r2_fq_path).unlink()

    with open(log_path, 'a') as log:
        log.write(f"\nAlignment complete: {output_bam}\n")
        log.write(f"Mode: {mode}\n")
        if use_r2:
            log.write(f"PAIRED-END: R1 and R2 aligned together\n")
        elif any_have_r2:
            log.write(f"MIXED: Some reads have R2, using R1 only for consistency\n")
        else:
            log.write(f"SINGLE-END: Only R1 aligned\n")


if __name__ == '__main__':
    main()
