"""
Command-line interface for TRACE.

TRACE: Triple-aligner Read Analysis for CRISPR Editing

Author: Kevin R. Roy
"""

import sys
from pathlib import Path

import click

from . import __version__
from .config import LocusConfig, NucleaseType, parse_sequence_input


@click.group()
@click.version_option(version=__version__)
def cli():
    """TRACE: Triple-aligner Read Analysis for CRISPR Editing."""
    pass


@cli.command()
@click.option('--reference', '-r', type=str,
              help='Reference amplicon: DNA sequence or FASTA file path (default for all samples)')
@click.option('--hdr-template', '-h', type=str,
              help='HDR template: DNA sequence or FASTA file path (default for all samples)')
@click.option('--guide', '-g', type=str,
              help='Guide sequence (20bp for Cas9, 20-24bp for Cas12a) (default for all samples)')
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

    REFERENCE, HDR_TEMPLATE, and GUIDE can be provided:
      - As CLI options (default for all samples)
      - Per-sample in the sample key TSV (reference, hdr_template, guide columns)

    Values can be DNA sequences directly or paths to FASTA files.

    \b
    Example with DNA sequences (single sample):
      trace run -r ATCG...250bp...ATCG -h ATCG...150bp...ATCG -g GCTGAAGCACTGCACGCCGT \\
                --r1 reads_R1.fastq.gz -o results/

    \b
    Example with FASTA files:
      trace run -r reference.fasta -h template.fasta -g GCTGAAGCACTGCACGCCGT \\
                --r1 reads_R1.fastq.gz -o results/

    \b
    Example with per-sample loci (sample key with reference/hdr_template/guide columns):
      trace run --sample-key samples.tsv -o results/
    """
    import logging

    from .config import PipelineConfig
    from .io.sample_key import Sample, load_sample_key
    from .pipeline import EditingPipeline
    from .preprocessing import run_auto_detection

    # Set up logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

    # Validate inputs
    if not r1 and not sample_key:
        click.echo("Error: Either --r1 or --sample-key must be provided", err=True)
        sys.exit(1)

    # For single sample mode, require sequences
    if r1 and not sample_key:
        if not all([reference, hdr_template, guide]):
            click.echo("Error: --reference, --hdr-template, and --guide are required for single sample mode", err=True)
            sys.exit(1)

    # Parse default sequences (if provided)
    default_ref_seq = None
    default_hdr_seq = None
    default_guide = guide

    if reference:
        try:
            default_ref_seq = parse_sequence_input(reference)
        except ValueError as e:
            click.echo(f"Error loading reference sequence: {e}", err=True)
            sys.exit(1)

    if hdr_template:
        try:
            default_hdr_seq = parse_sequence_input(hdr_template)
        except ValueError as e:
            click.echo(f"Error loading HDR template: {e}", err=True)
            sys.exit(1)

    nuclease_type = NucleaseType.CAS9 if nuclease == 'cas9' else NucleaseType.CAS12A

    # Load contaminant sequence if provided
    contaminant_seq = None
    if contaminant:
        try:
            contaminant_seq = parse_sequence_input(contaminant)
        except ValueError as e:
            click.echo(f"Error loading contaminant sequence: {e}", err=True)
            sys.exit(1)

    # Create output directory
    output_path = Path(output)
    output_path.mkdir(parents=True, exist_ok=True)

    # Handle single sample vs batch mode
    if r1:
        # Single sample mode - use CLI defaults
        ref_seq = default_ref_seq
        hdr_seq = default_hdr_seq

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

        locus.print_summary()

        # Create sample
        r1_path = Path(r1)
        r2_path = Path(r2) if r2 else None

        # Run auto-detection
        detection = run_auto_detection(r1_path, r2_path, len(ref_seq), reference=ref_seq)
        detection.print_summary()

        samples = [Sample(
            sample_id="sample",
            r1_path=r1_path,
            r2_path=r2_path,
        )]

        click.echo(f"\nProcessing single sample with {threads} threads...")

    else:
        # Batch mode - load samples
        samples = load_sample_key(Path(sample_key), validate=True)
        n_samples = len(samples)

        # Check if samples have per-sample loci or need CLI defaults
        samples_with_custom = [s for s in samples if s.has_custom_locus()]

        if samples_with_custom:
            # Per-sample loci detected - not yet fully supported in pipeline
            click.echo(f"\nNote: {len(samples_with_custom)} samples have per-sample locus overrides")
            click.echo("  Per-sample loci will be used for those samples.")

        # Need CLI defaults for samples without per-sample loci, or for the pipeline default
        if not all([default_ref_seq, default_hdr_seq, default_guide]):
            samples_without_custom = [s for s in samples if not s.has_custom_locus()]
            if samples_without_custom:
                missing_samples = [s.sample_id for s in samples_without_custom[:5]]
                click.echo(f"Error: {len(samples_without_custom)} samples need defaults", err=True)
                click.echo(f"  Provide --reference, --hdr-template, --guide for samples: {', '.join(missing_samples)}{'...' if len(samples_without_custom) > 5 else ''}", err=True)
                sys.exit(1)

        # Use CLI defaults for locus config
        ref_seq = default_ref_seq
        hdr_seq = default_hdr_seq

        try:
            locus = LocusConfig(
                name="analysis",
                reference=ref_seq,
                hdr_template=hdr_seq,
                guide=default_guide,
                nuclease=nuclease_type,
            ).analyze()
        except ValueError as e:
            click.echo(f"Error analyzing locus: {e}", err=True)
            sys.exit(1)

        locus.print_summary()

        click.echo(f"\nProcessing {n_samples} samples with {threads} threads...")

        # Run auto-detection on first sample
        if samples:
            first = samples[0]
            detection = run_auto_detection(first.r1_path, first.r2_path, len(ref_seq), reference=ref_seq)
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

    click.echo("\nPipeline complete!")
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


@cli.command('generate-manifest')
@click.option('--sample-key', '-s', type=click.Path(exists=True), required=True,
              help='Path to sample_key.tsv')
@click.option('--plate-key', '-p', type=click.Path(exists=True), required=True,
              help='Path to plate_key.tsv')
@click.option('--raw-data-dir', '-r', type=click.Path(exists=True), required=True,
              help='Path to raw_data directory containing FASTQ files')
@click.option('--output', '-o', type=click.Path(), required=True,
              help='Output path for TRACE sample manifest')
@click.option('--locus-filter', type=str, default=None,
              help='Filter to specific locus')
@click.option('--replicate-filter', type=str, default=None,
              help='Filter to specific replicate')
def generate_manifest(sample_key, plate_key, raw_data_dir, output, locus_filter, replicate_filter):
    """
    Generate TRACE-compatible sample manifest from project keyfiles.

    Joins sample_key with plate_key to resolve FASTQ paths and creates
    a sample manifest suitable for TRACE analysis.

    \\b
    Example:
      trace generate-manifest \\
        --sample-key keyfiles/sample_key.tsv \\
        --plate-key keyfiles/plate_key.tsv \\
        --raw-data-dir raw_data/ \\
        --output trace_sample_key.tsv
    """
    import logging

    from .utils.manifest import generate_trace_manifest

    logging.basicConfig(level=logging.INFO)

    try:
        result = generate_trace_manifest(
            sample_key_path=Path(sample_key),
            plate_key_path=Path(plate_key),
            raw_data_dir=Path(raw_data_dir),
            output_path=Path(output),
            locus_filter=locus_filter,
            replicate_filter=replicate_filter,
        )
        click.echo(f"Generated manifest with {len(result)} samples: {output}")
    except Exception as e:
        click.echo(f"Error generating manifest: {e}", err=True)
        sys.exit(1)


@cli.command('generate-templates')
@click.option('--sample-key', '-s', type=click.Path(exists=True),
              help='Path to sample_key.tsv (to extract barcodes)')
@click.option('--seq-ref', type=click.Path(exists=True),
              help='Path to guide_donor_and_reference_info.tsv (for homology arms)')
@click.option('--left-arm', type=str,
              help='Left homology arm sequence (alternative to --seq-ref)')
@click.option('--right-arm', type=str,
              help='Right homology arm sequence (alternative to --seq-ref)')
@click.option('--barcode-column', type=str, default='barcode',
              help='Column name containing barcodes in sample_key')
@click.option('--output', '-o', type=click.Path(), required=True,
              help='Output FASTA path for HDR templates')
def generate_templates(sample_key, seq_ref, left_arm, right_arm, barcode_column, output):
    """
    Generate HDR templates FASTA with all possible barcode insertions.

    Can extract homology arms from guide_donor_and_reference_info.tsv or
    accept them directly via --left-arm and --right-arm.

    \\b
    Example using keyfiles:
      trace generate-templates \\
        --sample-key keyfiles/sample_key.tsv \\
        --seq-ref keyfiles/guide_donor_and_reference_info.tsv \\
        --output hdr_templates.fasta

    \\b
    Example with direct arms:
      trace generate-templates \\
        --sample-key keyfiles/sample_key.tsv \\
        --left-arm "ATCGATCG..." \\
        --right-arm "GCTAGCTA..." \\
        --output hdr_templates.fasta
    """
    import logging

    from .utils.barcode_templates import (
        generate_barcoded_hdr_templates,
        generate_templates_from_keyfiles,
        load_barcodes_from_tsv,
    )

    logging.basicConfig(level=logging.INFO)

    try:
        if seq_ref:
            # Use keyfiles for homology arms
            if not sample_key:
                click.echo("Error: --sample-key required when using --seq-ref", err=True)
                sys.exit(1)

            templates = generate_templates_from_keyfiles(
                sample_key_path=Path(sample_key),
                guide_donor_info_path=Path(seq_ref),
                output_fasta=Path(output),
                barcode_column=barcode_column,
            )
        elif left_arm and right_arm:
            # Use direct arm sequences
            if not sample_key:
                click.echo("Error: --sample-key required to extract barcodes", err=True)
                sys.exit(1)

            barcodes = load_barcodes_from_tsv(Path(sample_key), barcode_column)
            templates = generate_barcoded_hdr_templates(
                left_homology_arm=left_arm,
                right_homology_arm=right_arm,
                barcodes=barcodes,
                output_fasta=Path(output),
            )
        else:
            click.echo("Error: Provide either --seq-ref or both --left-arm and --right-arm", err=True)
            sys.exit(1)

        click.echo(f"Generated {len(templates)} HDR templates: {output}")

    except Exception as e:
        click.echo(f"Error generating templates: {e}", err=True)
        import traceback
        traceback.print_exc()
        sys.exit(1)


@cli.command('multi-template')
@click.option('--reference', '-r', type=str, required=True,
              help='Reference amplicon: DNA sequence or FASTA file path')
@click.option('--hdr-templates', '-H', type=click.Path(exists=True), required=True,
              help='FASTA file with multiple HDR templates (one per barcode)')
@click.option('--guide', '-g', type=str, required=True,
              help='Guide sequence')
@click.option('--sample-key', '-s', type=click.Path(exists=True), required=True,
              help='Sample key TSV file')
@click.option('--output', '-o', type=click.Path(), required=True,
              help='Output directory')
@click.option('--nuclease', type=click.Choice(['cas9', 'cas12a']), default='cas9',
              help='Nuclease type (default: cas9)')
@click.option('--expected-template-column', type=str, default='expected_barcode',
              help='Column in sample key with expected template for each sample')
@click.option('--threads', '-t', type=int, default=4,
              help='Number of threads (default: 4)')
def multi_template(reference, hdr_templates, guide, sample_key, output,
                   nuclease, expected_template_column, threads):
    """
    Run multi-template analysis for barcode screening experiments.

    Analyzes samples for the presence of multiple possible HDR templates,
    allowing detection of expected and unexpected barcodes.

    \\b
    Example:
      trace multi-template \\
        --reference reference.fasta \\
        --hdr-templates hdr_templates.fasta \\
        --guide GAGTCCGAGCAGAAGAAGAA \\
        --sample-key trace_sample_key.tsv \\
        --expected-template-column expected_barcode \\
        --output results/
    """
    import logging

    from .config import MultiTemplateLocusConfig, NucleaseType
    from .io.sample_key import load_sample_key
    from .pipeline import run_multi_template_pipeline

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

    # Load reference sequence
    try:
        ref_seq = parse_sequence_input(reference)
    except ValueError as e:
        click.echo(f"Error loading reference: {e}", err=True)
        sys.exit(1)

    # Load HDR templates from FASTA
    try:
        templates = {}
        current_id = None
        current_seq = []
        with open(hdr_templates) as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_id:
                        templates[current_id] = ''.join(current_seq).upper()
                    current_id = line[1:].split()[0]
                    current_seq = []
                else:
                    current_seq.append(line)
            if current_id:
                templates[current_id] = ''.join(current_seq).upper()

        click.echo(f"Loaded {len(templates)} HDR templates")

    except Exception as e:
        click.echo(f"Error loading HDR templates: {e}", err=True)
        sys.exit(1)

    # Load samples
    try:
        samples = load_sample_key(Path(sample_key), validate=True)
        click.echo(f"Loaded {len(samples)} samples")
    except Exception as e:
        click.echo(f"Error loading sample key: {e}", err=True)
        sys.exit(1)

    nuclease_type = NucleaseType.CAS9 if nuclease == 'cas9' else NucleaseType.CAS12A

    # Create and analyze locus
    try:
        locus = MultiTemplateLocusConfig(
            name="analysis",
            reference=ref_seq,
            hdr_templates=templates,
            guide=guide,
            nuclease=nuclease_type,
        ).analyze()

        locus.print_summary()

    except ValueError as e:
        click.echo(f"Error analyzing locus: {e}", err=True)
        sys.exit(1)

    # Create output directory
    output_path = Path(output)
    output_path.mkdir(parents=True, exist_ok=True)

    # Run pipeline
    click.echo(f"\nRunning multi-template pipeline with {threads} threads...")

    results = run_multi_template_pipeline(
        reference=ref_seq,
        hdr_templates=templates,
        guide=guide,
        samples=samples,
        output_dir=output_path,
        nuclease=nuclease,
        expected_template_column=expected_template_column,
        threads=threads,
    )

    click.echo("\nPipeline complete!")
    click.echo(f"Processed {len(results)} samples")
    click.echo(f"Results written to: {output_path}")


if __name__ == '__main__':
    cli()
