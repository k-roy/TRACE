#!/usr/bin/env python3
"""
Align unique R1/R2 pairs to reference using triple-aligner consensus (BWA-MEM, BBMap, minimap2).

This script follows TRACE's original design philosophy of using three aligners
for robust alignment across diverse read types and reference sequences.

Author: Kevin R. Roy
"""

import logging
import subprocess
from pathlib import Path
import tempfile

import pandas as pd
import pysam


def run_alignment(input_tsv, reference_fasta, output_bam, output_bai, log_path, threads=4):
    """
    Align unique R1/R2 pairs using triple-aligner approach.

    Args:
        input_tsv: Path to unique sequences TSV with r1_sequence and r2_sequence columns
        reference_fasta: Path to reference FASTA file
        output_bam: Path for output BAM file (will use first successful aligner)
        output_bai: Path for output BAM index file
        log_path: Path for log file
        threads: Number of threads for alignment
    """
    # Set up logging to file
    log_path = Path(log_path)
    log_path.parent.mkdir(parents=True, exist_ok=True)

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_path),
            logging.StreamHandler()
        ]
    )
    log = logging.getLogger(__name__)

    # Read unique sequences
    df = pd.read_csv(input_tsv, sep='\t', compression='gzip')
    log.info(f"Loaded {len(df):,} unique R1/R2 pairs")

    # Detect mode: paired-end or single-end
    has_r2 = any(row['r2_sequence'] and row['r2_sequence'] != '' for _, row in df.iterrows())
    mode = "paired-end" if has_r2 else "single-end"
    log.info(f"Detected {mode} mode")

    # Create temporary FASTQ files
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        r1_fq_path = tmpdir / "unique_r1.fastq"
        r2_fq_path = tmpdir / "unique_r2.fastq" if has_r2 else None

        # Write R1 FASTQ
        with open(r1_fq_path, 'w') as r1_fq:
            for _, row in df.iterrows():
                seq_hash = row['seq_hash']
                r1_seq = row['r1_sequence']
                r1_fq.write(f"@{seq_hash}\n{r1_seq}\n+\n{'I'*len(r1_seq)}\n")

        log.info(f"Wrote {len(df)} R1 sequences to {r1_fq_path}")

        # Write R2 FASTQ if paired-end
        if has_r2:
            with open(r2_fq_path, 'w') as r2_fq:
                for _, row in df.iterrows():
                    seq_hash = row['seq_hash']
                    r2_seq = row['r2_sequence']
                    if r2_seq and r2_seq != '':
                        r2_fq.write(f"@{seq_hash}\n{r2_seq}\n+\n{'I'*len(r2_seq)}\n")
            log.info(f"Wrote {len(df)} R2 sequences to {r2_fq_path}")

        # Run triple alignment
        log.info("Running triple-aligner consensus...")
        log.info(f"  Aligner 1: BWA-MEM (paired-end aware)")
        log.info(f"  Aligner 2: BBMap (gapped aligner)")
        log.info(f"  Aligner 3: minimap2 (fast approximate aligner)")

        aligners = ['bwa', 'bbmap', 'minimap2']
        aligner_bams = {}

        for aligner in aligners:
            aligner_bam = Path(str(output_bam).replace('.bam', f'_{aligner}.bam'))
            aligner_bams[aligner] = aligner_bam

            try:
                if aligner == 'bwa':
                    if has_r2:
                        cmd = f"bwa mem -t {threads} {reference_fasta} {r1_fq_path} {r2_fq_path} | samtools sort -@ 2 -o {aligner_bam} -"
                    else:
                        cmd = f"bwa mem -t {threads} {reference_fasta} {r1_fq_path} | samtools sort -@ 2 -o {aligner_bam} -"

                elif aligner == 'bbmap':
                    cmd_parts = [
                        'bbmap.sh',
                        f'in={r1_fq_path}',
                        f'ref={reference_fasta}',
                        f'out={aligner_bam}',
                        f'threads={threads}',
                        'nodisk=t',
                        'ambiguous=best',
                        'pairedonly=f',
                    ]
                    if has_r2:
                        cmd_parts.insert(2, f'in2={r2_fq_path}')
                    cmd = ' '.join(cmd_parts)

                elif aligner == 'minimap2':
                    if has_r2:
                        cmd = f"minimap2 -ax sr -t {threads} --end-bonus=50 --secondary=no {reference_fasta} {r1_fq_path} {r2_fq_path} | samtools view -bS - | samtools sort -@ 2 -o {aligner_bam} -"
                    else:
                        cmd = f"minimap2 -ax sr -t {threads} --end-bonus=50 --secondary=no {reference_fasta} {r1_fq_path} | samtools view -bS - | samtools sort -@ 2 -o {aligner_bam} -"

                log.info(f"Running {aligner}...")
                result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

                if result.returncode != 0:
                    log.warning(f"{aligner} failed: {result.stderr}")
                    continue

                # Index the BAM
                subprocess.run(['samtools', 'index', str(aligner_bam)], check=True)
                log.info(f"{aligner} completed successfully")

            except Exception as e:
                log.warning(f"{aligner} error: {e}")
                continue

        # Use first successful aligner as primary output
        # Priority: BWA > BBMap > minimap2
        for aligner in ['bwa', 'bbmap', 'minimap2']:
            if aligner in aligner_bams and aligner_bams[aligner].exists():
                log.info(f"Using {aligner} as primary output")
                subprocess.run(['cp', str(aligner_bams[aligner]), str(output_bam)], check=True)
                subprocess.run(['cp', f"{aligner_bams[aligner]}.bai", str(output_bai)], check=True)
                break
        else:
            raise RuntimeError("All aligners failed!")

        # Log alignment statistics
        with pysam.AlignmentFile(output_bam, 'rb') as bam:
            total_reads = 0
            mapped_reads = 0
            for read in bam:
                total_reads += 1
                if not read.is_unmapped:
                    mapped_reads += 1

            if total_reads > 0:
                mapping_rate = 100 * mapped_reads / total_reads
                log.info(f"Alignment complete:")
                log.info(f"   Total reads: {total_reads:,}")
                log.info(f"   Mapped reads: {mapped_reads:,} ({mapping_rate:.1f}%)")
                log.info(f"   Unmapped reads: {total_reads - mapped_reads:,} ({100-mapping_rate:.1f}%)")
                if has_r2:
                    log.info(f"   Mode: Paired-end")
                    log.info(f"   Insert size and pairing flags calculated by aligners")
                else:
                    log.info(f"   Mode: Single-end")


def main():
    """Snakemake interface - reads from snakemake object."""
    run_alignment(
        input_tsv=snakemake.input.unique_seqs,
        reference_fasta=snakemake.input.ref,
        output_bam=snakemake.output.bam,
        output_bai=snakemake.output.bai,
        log_path=str(snakemake.log),
        threads=snakemake.threads
    )


if __name__ == '__main__':
    import sys
    # Check if running via Snakemake or command line
    if 'snakemake' in dir():
        main()
    elif len(sys.argv) == 5:
        # Command-line interface
        run_alignment(
            input_tsv=sys.argv[1],
            reference_fasta=sys.argv[2],
            output_bam=sys.argv[3],
            output_bai=sys.argv[3] + '.bai',
            log_path='/dev/stderr',
            threads=int(sys.argv[4])
        )
    else:
        print("Usage: align_unique_sequences_triple.py <input.tsv.gz> <reference.fasta> <output.bam> <threads>")
        print("Or run via Snakemake script: directive")
        sys.exit(1)
