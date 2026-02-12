"""
Command-line interface for TRACE.

TRACE: Triple-aligner Read Analysis for CRISPR Editing

Author: Kevin R. Roy
"""

import click
from pathlib import Path
from typing import Optional
import sys

from . import __version__
from .config import LocusConfig, NucleaseType, parse_sequence_input


@click.group()
@click.version_option(version=__version__)
def cli():
    """TRACE: Triple-aligner Read Analysis for CRISPR Editing."""
    pass


@cli.command()
@click.option('--reference', '-r', type=str, required=True,
              help='Reference amplicon: DNA sequence or FASTA file path')
@click.option('--hdr-template', '-h', type=str, required=True,
              help='HDR template: DNA sequence or FASTA file path (can be shorter than reference)')
@click.option('--guide', '-g', type=str, required=True,
              help='Guide sequence (20bp for Cas9, 20-24bp for Cas12a)')
@click.option('--r1', type=click.Path(exists=True),
              help='R1 FASTQ file (for single sample)')
@click.option('--r2', type=click.Path(exists=True),
              help='R2 FASTQ file (for paired-end)')
@click.option('--sample-key', '-s', type=click.Path(exists=True),
              help='Sample key TSV file (for multiple samples)')
@click.option('--output', '-o', type=click.Path(), required=True,
              help='Output directory')
@click.option('--nuclease', type=click.Choice(['cas9', 'cas12a']), default='cas9',
              help='Nuclease type (default: cas9)')
@click.option('--contaminant', '-c', type=str,
              help='Optional contaminant: DNA sequence or FASTA file path')
@click.option('--threads', '-t', type=int, default=4,
              help='Number of threads (default: 4)')
@click.option('--crispresso/--no-crispresso', default=True,
              help='Run CRISPResso2 for comparison (default: enabled)')
def run(reference, hdr_template, guide, r1, r2, sample_key, output,
        nuclease, contaminant, threads, crispresso):
    """
    Run the full editing outcome analysis pipeline.

    REFERENCE and HDR_TEMPLATE can be provided as:
      - DNA sequences directly (e.g., ATCGATCG...)
      - Paths to FASTA files (e.g., reference.fasta)

    \b
    Example with DNA sequences:
      trace run -r ATCG...250bp...ATCG -h ATCG...150bp...ATCG -g GCTGAAGCACTGCACGCCGT \\
                --r1 reads_R1.fastq.gz -o results/

    \b
    Example with FASTA files:
      trace run -r reference.fasta -h template.fasta -g GCTGAAGCACTGCACGCCGT \\
                --r1 reads_R1.fastq.gz -o results/
    """

    # Validate inputs
    if not r1 and not sample_key:
        click.echo("Error: Either --r1 or --sample-key must be provided", err=True)
        sys.exit(1)

    # Parse sequences (handles both DNA strings and file paths)
    try:
        ref_seq = parse_sequence_input(reference)
        hdr_seq = parse_sequence_input(hdr_template)
    except ValueError as e:
        click.echo(f"Error loading sequences: {e}", err=True)
        sys.exit(1)

    # Create locus configuration
    nuclease_type = NucleaseType.CAS9 if nuclease == 'cas9' else NucleaseType.CAS12A

    try:
        locus = LocusConfig(
            name="analysis",
            reference=ref_seq,
            hdr_template=hdr_seq,
            guide=guide,
            nuclease=nuclease_type,
        ).analyze()
    except ValueError as e:
        click.echo(f"Error analyzing locus: {e}", err=True)
        sys.exit(1)

    # Print configuration summary
    locus.print_summary()

    # Create output directory
    output_path = Path(output)
    output_path.mkdir(parents=True, exist_ok=True)

    # Import pipeline components
    from .preprocessing import run_auto_detection
    from .io.sample_key import Sample, load_sample_key
    from .pipeline import EditingPipeline
    from .config import PipelineConfig
    import logging

    # Set up logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

    # Load contaminant sequence if provided
    contaminant_seq = None
    if contaminant:
        try:
            contaminant_seq = parse_sequence_input(contaminant)
        except ValueError as e:
            click.echo(f"Error loading contaminant sequence: {e}", err=True)
            sys.exit(1)

    # Handle single sample vs batch mode
    if r1:
        # Single sample mode
        r1_path = Path(r1)
        r2_path = Path(r2) if r2 else None

        # Run auto-detection
        detection = run_auto_detection(r1_path, r2_path, len(ref_seq))
        detection.print_summary()

        # Create sample
        samples = [Sample(
            sample_id="sample",
            r1_path=r1_path,
            r2_path=r2_path,
        )]

        click.echo(f"\nProcessing single sample with {threads} threads...")

    else:
        # Batch mode with sample key
        samples = load_sample_key(Path(sample_key), validate=True)
        n_samples = len(samples)

        click.echo(f"\nProcessing {n_samples} samples with {threads} threads...")

        # Get first sample for auto-detection preview
        if samples:
            first = samples[0]
            detection = run_auto_detection(first.r1_path, first.r2_path, len(ref_seq))
            detection.print_summary()

    # Create pipeline configuration
    config = PipelineConfig(
        locus=locus,
        samples=samples,
        output_dir=output_path,
        threads=threads,
        run_crispresso=crispresso,
        contaminant_sequence=contaminant_seq,
    )

    # Run pipeline
    click.echo("\nRunning pipeline...")
    pipeline = EditingPipeline(config)
    results = pipeline.run()

    click.echo(f"\nPipeline complete!")
    click.echo(f"Processed {len(results)} samples")
    click.echo(f"Results written to: {output_path / 'per_sample_editing_outcomes_all_methods.tsv'}")


@cli.command()
@click.argument('fastq', type=click.Path(exists=True), nargs=-1, required=True)
@click.option('--reference', '-r', type=str, required=True,
              help='Reference amplicon: DNA sequence or FASTA file path')
@click.option('--output', '-o', type=click.Path(), required=True,
              help='Output BAM file')
@click.option('--threads', '-t', type=int, default=4,
              help='Number of threads')
def align(fastq, reference, output, threads):
    """Align reads using triple-aligner approach (BWA, BBMap, minimap2)."""
    click.echo(f"Aligning {len(fastq)} FASTQ file(s)")
    click.echo(f"Output: {output}")
    click.echo(f"Threads: {threads}")

    # TODO: Implement alignment
    click.echo("\nAlignment not yet implemented.")


@cli.command()
@click.argument('bam', type=click.Path(exists=True))
@click.option('--reference', '-r', type=str, required=True,
              help='Reference amplicon: DNA sequence or FASTA file path')
@click.option('--hdr-template', '-h', type=str, required=True,
              help='HDR template: DNA sequence or FASTA file path')
@click.option('--guide', '-g', type=str, required=True,
              help='Guide sequence')
@click.option('--output', '-o', type=click.Path(), required=True,
              help='Output TSV file')
@click.option('--nuclease', type=click.Choice(['cas9', 'cas12a']), default='cas9',
              help='Nuclease type')
def classify(bam, reference, hdr_template, guide, output, nuclease):
    """Classify aligned reads into editing outcomes."""
    click.echo(f"Classifying reads from: {bam}")

    # Load sequences
    try:
        ref_seq = parse_sequence_input(reference)
        hdr_seq = parse_sequence_input(hdr_template)
    except ValueError as e:
        click.echo(f"Error loading sequences: {e}", err=True)
        sys.exit(1)

    nuclease_type = NucleaseType.CAS9 if nuclease == 'cas9' else NucleaseType.CAS12A

    locus = LocusConfig(
        name="analysis",
        reference=ref_seq,
        hdr_template=hdr_seq,
        guide=guide,
        nuclease=nuclease_type,
    ).analyze()

    locus.print_summary()

    # TODO: Implement classification
    click.echo("\nClassification not yet implemented.")


@cli.command('extract-kmers')
@click.argument('contaminant', type=str)
@click.option('--reference', '-r', type=str,
              help='Reference sequence to exclude from k-mers (DNA or FASTA path)')
@click.option('--kmer-size', '-k', type=int, default=12,
              help='K-mer size (default: 12)')
@click.option('--output', '-o', type=click.Path(), required=True,
              help='Output file for k-mers')
def extract_kmers(contaminant, reference, kmer_size, output):
    """Extract unique k-mers from contaminant sequence for filtering."""
    from .utils.sequence import extract_unique_kmers

    try:
        contaminant_seq = parse_sequence_input(contaminant)
    except ValueError as e:
        click.echo(f"Error loading contaminant sequence: {e}", err=True)
        sys.exit(1)

    exclude_seqs = []
    if reference:
        try:
            exclude_seqs.append(parse_sequence_input(reference))
        except ValueError as e:
            click.echo(f"Error loading reference sequence: {e}", err=True)
            sys.exit(1)

    kmers = extract_unique_kmers(contaminant_seq, kmer_size, exclude_seqs)

    # Write k-mers to file
    with open(output, 'w') as f:
        for kmer in sorted(kmers):
            f.write(f"{kmer}\n")

    click.echo(f"Extracted {len(kmers)} unique {kmer_size}-mers to {output}")


@cli.command()
@click.option('--output', '-o', type=click.Path(), default='trace_config.yaml',
              help='Output config file path')
def init(output):
    """Generate a template configuration file."""
    template = '''# TRACE Configuration Template
# Edit this file to configure your analysis

# Required: Sequence files
reference: amplicon.fasta           # Reference amplicon sequence
hdr_template: hdr_template.fasta    # HDR template sequence
guide: GCTGAAGCACTGCACGCCGT         # Guide sequence (without PAM)

# Nuclease type: cas9 or cas12a
nuclease: cas9

# Sample information (use sample_key for multiple samples)
# For single sample:
# r1: sample_R1.fastq.gz
# r2: sample_R2.fastq.gz

# For multiple samples, create a TSV file:
sample_key: samples.tsv

# Output directory
output_dir: ./results

# Processing options
threads: 16

# Optional: Contaminant filtering
# contaminant: gfp_plasmid.fasta

# CRISPResso2 integration (enabled by default)
crispresso: true
'''

    with open(output, 'w') as f:
        f.write(template)

    click.echo(f"Generated configuration template: {output}")
    click.echo("\nEdit this file and run:")
    click.echo(f"  trace run --config {output}")


@cli.command()
@click.option('--reference', '-r', type=str, required=True,
              help='Reference amplicon: DNA sequence or FASTA file path')
@click.option('--hdr-template', '-h', type=str, required=True,
              help='HDR template: DNA sequence or FASTA file path')
@click.option('--guide', '-g', type=str, required=True,
              help='Guide sequence')
@click.option('--nuclease', type=click.Choice(['cas9', 'cas12a']), default='cas9',
              help='Nuclease type')
def info(reference, hdr_template, guide, nuclease):
    """
    Display locus analysis information without running the pipeline.

    \b
    Example with DNA sequences (250bp reference, 150bp template):
      trace info \\
        -r ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG... \\
        -h ATCGATCGATCGATCGATCGATCGATCGATCGATCG... \\
        -g GCTGAAGCACTGCACGCCGT

    \b
    Example with FASTA files:
      trace info -r reference.fasta -h template.fasta -g GCTGAAGCACTGCACGCCGT
    """
    try:
        ref_seq = parse_sequence_input(reference)
        hdr_seq = parse_sequence_input(hdr_template)
    except ValueError as e:
        click.echo(f"Error loading sequences: {e}", err=True)
        sys.exit(1)

    nuclease_type = NucleaseType.CAS9 if nuclease == 'cas9' else NucleaseType.CAS12A

    try:
        locus = LocusConfig(
            name="analysis",
            reference=ref_seq,
            hdr_template=hdr_seq,
            guide=guide,
            nuclease=nuclease_type,
        ).analyze()

        locus.print_summary()

    except ValueError as e:
        click.echo(f"Error analyzing locus: {e}", err=True)
        sys.exit(1)


if __name__ == '__main__':
    cli()
