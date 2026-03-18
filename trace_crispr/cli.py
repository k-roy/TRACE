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
        samples_without_custom = [s for s in samples if not s.has_custom_locus()]

        if samples_with_custom:
            click.echo(f"\n✓ {len(samples_with_custom)}/{n_samples} samples have per-sample sequences")

            # Validate that per-sample sequences are complete
            incomplete = []
            for s in samples_with_custom:
                missing_fields = []
                if not s.reference and not default_ref_seq:
                    missing_fields.append('reference')
                if not s.guide and not default_guide:
                    missing_fields.append('guide')
                if not s.hdr_template and not default_hdr_seq:
                    missing_fields.append('hdr_template')

                if missing_fields:
                    incomplete.append(f"{s.sample_id}: missing {', '.join(missing_fields)}")

            if incomplete:
                click.echo("\nERROR: Samples with incomplete per-sample sequences:", err=True)
                for msg in incomplete[:10]:  # Show first 10
                    click.echo(f"  - {msg}", err=True)
                if len(incomplete) > 10:
                    click.echo(f"  ... and {len(incomplete)-10} more", err=True)
                sys.exit(1)
        else:
            click.echo(f"\nUsing global sequences for all {n_samples} samples")

        # Check if defaults are needed
        if samples_without_custom and not all([default_ref_seq, default_hdr_seq, default_guide]):
            missing_samples = [s.sample_id for s in samples_without_custom[:5]]
            click.echo(
                f"\nERROR: {len(samples_without_custom)} samples don't have per-sample sequences, "
                "but no defaults provided via --reference, --guide, --hdr-template",
                err=True
            )
            click.echo(f"  Samples needing defaults: {', '.join(missing_samples)}{'...' if len(samples_without_custom) > 5 else ''}", err=True)
            sys.exit(1)

        # Use CLI defaults for locus config (may be None if all samples have custom sequences)
        ref_seq = default_ref_seq
        hdr_seq = default_hdr_seq

        # Create default locus config (if defaults provided)
        if ref_seq and hdr_seq and default_guide:
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
        else:
            # All samples must have per-sample sequences
            # Create a dummy locus using the first sample's sequences for initialization
            first_custom_sample = samples_with_custom[0]
            locus = LocusConfig(
                name="default",
                reference=first_custom_sample.reference or "ATCG",  # Dummy if needed
                hdr_template=first_custom_sample.hdr_template or "ATCG",
                guide=first_custom_sample.guide or "ATCGATCGATCGATCGATCG",
                nuclease=nuclease_type,
            ).analyze()
            click.echo("\nNo global defaults - all samples use per-sample sequences")

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
@click.option('--per-read', is_flag=True, default=False,
              help='Also write per-read classification details')
def classify(bam, reference, hdr_template, guide, output, nuclease, per_read):
    """
    Classify aligned reads into editing outcomes.

    Re-run classification from existing BAM files without re-running alignment.
    Useful when classification scheme has been updated.

    \b
    Example:
      trace classify alignments/sample.bam \\
        -r reference.fasta \\
        -h hdr_template.fasta \\
        -g GCTGAAGCACTGCACGCCGT \\
        -o classified_results.tsv
    """
    import logging

    import pysam

    from .core.classification import (
        EditingOutcome,
        classify_read,
        get_hdr_signature_positions,
        summarize_classifications,
    )
    from .io.output import write_per_read_classifications

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

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

    # Get HDR signature
    min_len = min(len(ref_seq), len(hdr_seq))
    hdr_signature = get_hdr_signature_positions(ref_seq[:min_len], hdr_seq[:min_len])
    click.echo(f"HDR signature: {len(hdr_signature)} positions differ from WT")

    # Get cut site
    cut_site = locus.guide_info.cleavage_site if locus.guide_info else 0

    # Classify reads
    classifications = []
    per_read_details = []
    aligned_reads = 0
    unmapped_reads = 0

    click.echo("Classifying reads...")
    try:
        with pysam.AlignmentFile(bam, "rb") as bam_file:
            for read in bam_file.fetch(until_eof=True):
                if read.is_unmapped:
                    unmapped_reads += 1
                    continue

                aligned_reads += 1

                result = classify_read(
                    read,
                    hdr_signature,
                    cut_site,
                    ref_seq=ref_seq,
                    donor_seq=hdr_seq,
                )

                classifications.append(result)

                if per_read:
                    per_read_details.append({
                        'read_name': read.query_name,
                        'outcome': result.outcome.value,
                        'confidence': f"{result.confidence:.3f}",
                        'qc_flags': ','.join(result.qc_flags) if result.qc_flags else '',
                        **{k: str(v) for k, v in result.details.items() if not isinstance(v, (dict, list))}
                    })

    except Exception as e:
        click.echo(f"Error reading BAM: {e}", err=True)
        sys.exit(1)

    # Summarize
    summary = summarize_classifications(classifications)

    click.echo(f"\n--- Classification Summary ---")
    click.echo(f"Total aligned reads: {aligned_reads}")
    click.echo(f"Unmapped reads: {unmapped_reads}")
    click.echo(f"\n--- Outcome Breakdown ---")

    # Print each category
    for outcome in EditingOutcome:
        count = summary.get(f'{outcome.value}_count', 0)
        rate = summary.get(f'{outcome.value}_rate', 0) * 100
        if count > 0:
            click.echo(f"  {outcome.value}: {count} ({rate:.2f}%)")

    # Print aggregated metrics
    click.echo(f"\n--- Aggregated Metrics ---")
    click.echo(f"  HDR total: {summary.get('hdr_total_count', 0)} ({summary.get('hdr_total_rate', 0)*100:.2f}%)")
    click.echo(f"  NHEJ/MMEJ total: {summary.get('nhej_mmej_total_count', 0)} ({summary.get('nhej_mmej_total_rate', 0)*100:.2f}%)")
    click.echo(f"  Edited total: {summary.get('edited_total_count', 0)} ({summary.get('edited_total_rate', 0)*100:.2f}%)")

    # Write summary output
    output_path = Path(output)
    import pandas as pd

    summary_df = pd.DataFrame([{
        'bam_file': bam,
        'total_aligned': aligned_reads,
        'total_unmapped': unmapped_reads,
        **{f'{outcome.value}_count': summary.get(f'{outcome.value}_count', 0)
           for outcome in EditingOutcome},
        **{f'{outcome.value}_pct': summary.get(f'{outcome.value}_rate', 0) * 100
           for outcome in EditingOutcome},
        'hdr_total_count': summary.get('hdr_total_count', 0),
        'hdr_total_pct': summary.get('hdr_total_rate', 0) * 100,
        'nhej_mmej_total_count': summary.get('nhej_mmej_total_count', 0),
        'nhej_mmej_total_pct': summary.get('nhej_mmej_total_rate', 0) * 100,
        'edited_total_count': summary.get('edited_total_count', 0),
        'edited_total_pct': summary.get('edited_total_rate', 0) * 100,
    }])

    summary_df.to_csv(output_path, sep='\t', index=False)
    click.echo(f"\nSummary written to: {output_path}")

    # Write per-read details if requested
    if per_read and per_read_details:
        per_read_path = output_path.parent / f"{output_path.stem}_per_read.tsv"
        write_per_read_classifications(per_read_details, per_read_path)
        click.echo(f"Per-read details written to: {per_read_path}")


@cli.command('classify-batch')
@click.option('--bam-dir', '-b', type=click.Path(exists=True), required=True,
              help='Directory containing BAM files from previous TRACE run')
@click.option('--sample-key', '-s', type=click.Path(exists=True),
              help='Optional sample key TSV for metadata (uses BAM names if not provided)')
@click.option('--reference', '-r', type=str, required=True,
              help='Reference amplicon: DNA sequence or FASTA file path')
@click.option('--hdr-template', '-h', type=str, required=True,
              help='HDR template: DNA sequence or FASTA file path')
@click.option('--guide', '-g', type=str, required=True,
              help='Guide sequence')
@click.option('--output', '-o', type=click.Path(), required=True,
              help='Output directory')
@click.option('--nuclease', type=click.Choice(['cas9', 'cas12a']), default='cas9',
              help='Nuclease type')
@click.option('--bam-pattern', type=str, default='*_sorted.bam',
              help='Glob pattern for BAM files (default: *_sorted.bam)')
@click.option('--threads', '-t', type=int, default=4,
              help='Number of threads')
def classify_batch(bam_dir, sample_key, reference, hdr_template, guide, output,
                   nuclease, bam_pattern, threads):
    """
    Re-run classification on multiple samples from existing BAM files.

    Use this to re-classify reads with an updated classification scheme
    without re-running alignment. Perfect for iterating on classification
    logic after alignment is complete.

    \b
    Example:
      trace classify-batch \\
        --bam-dir results/sample1/alignments \\
        --sample-key samples.tsv \\
        -r reference.fasta \\
        -h hdr_template.fasta \\
        -g GCTGAAGCACTGCACGCCGT \\
        -o reclassified_results/
    """
    import logging
    from concurrent.futures import ProcessPoolExecutor, as_completed

    import pysam

    from .core.classification import (
        EditingOutcome,
        classify_read,
        get_hdr_signature_positions,
        summarize_classifications,
    )
    from .io.output import SampleResult, write_results_tsv, generate_summary_report
    from .io.sample_key import load_sample_key

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

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

    # Get HDR signature and cut site
    min_len = min(len(ref_seq), len(hdr_seq))
    hdr_signature = get_hdr_signature_positions(ref_seq[:min_len], hdr_seq[:min_len])
    cut_site = locus.guide_info.cleavage_site if locus.guide_info else 0

    click.echo(f"HDR signature: {len(hdr_signature)} positions")
    click.echo(f"Cut site: {cut_site}")

    # Find BAM files
    bam_dir_path = Path(bam_dir)
    bam_files = list(bam_dir_path.glob(bam_pattern))

    # Also check subdirectories (for per-sample output structure)
    if not bam_files:
        bam_files = list(bam_dir_path.glob(f"**/alignments/{bam_pattern}"))

    if not bam_files:
        click.echo(f"No BAM files found matching '{bam_pattern}' in {bam_dir}", err=True)
        sys.exit(1)

    click.echo(f"Found {len(bam_files)} BAM files")

    # Load sample metadata if provided
    sample_metadata = {}
    if sample_key:
        samples = load_sample_key(Path(sample_key), validate=False)
        sample_metadata = {s.sample_id: s.metadata for s in samples}
        click.echo(f"Loaded metadata for {len(samples)} samples")

    # Create output directory
    output_path = Path(output)
    output_path.mkdir(parents=True, exist_ok=True)

    # Process each BAM file
    results = []
    for i, bam_path in enumerate(bam_files):
        # Extract sample ID from BAM filename
        sample_id = bam_path.stem.replace('_sorted', '').replace('_bbmap', '').replace('_bwa', '')

        click.echo(f"Processing {i+1}/{len(bam_files)}: {sample_id}")

        # Classify reads
        classifications = []
        aligned_reads = 0
        unmapped_reads = 0

        try:
            with pysam.AlignmentFile(str(bam_path), "rb") as bam_file:
                for read in bam_file.fetch(until_eof=True):
                    if read.is_unmapped:
                        unmapped_reads += 1
                        continue

                    aligned_reads += 1

                    result = classify_read(
                        read,
                        hdr_signature,
                        cut_site,
                        ref_seq=ref_seq,
                        donor_seq=hdr_seq,
                    )
                    classifications.append(result)

        except Exception as e:
            click.echo(f"  Error reading BAM: {e}", err=True)
            continue

        # Summarize
        summary = summarize_classifications(classifications)

        # Build SampleResult
        sample_result = SampleResult(
            sample_id=sample_id,
            total_reads=aligned_reads + unmapped_reads,
            aligned_reads=aligned_reads,
            classifiable_reads=len(classifications),

            # New comprehensive classification counts
            dedup_hdr_complete_count=summary.get('hdr_complete_count', 0),
            dedup_hdr_partial_count=summary.get('hdr_partial_count', 0),
            dedup_hdr_plus_nhej_count=summary.get('hdr_plus_nhej_indel_count', 0),
            dedup_hdr_plus_mmej_count=summary.get('hdr_plus_mmej_indel_count', 0),
            dedup_hdr_plus_other_count=summary.get('hdr_plus_other_count', 0),
            dedup_donor_capture_count=summary.get('donor_capture_count', 0),
            dedup_nhej_indel_count=summary.get('nhej_indel_count', 0),
            dedup_mmej_indel_count=summary.get('mmej_indel_count', 0),
            dedup_wt_count=summary.get('wt_count', 0),
            dedup_non_donor_snv_count=summary.get('non_donor_snv_count', 0),
            dedup_unclassified_count=summary.get('unclassified_count', 0),
            dedup_unmapped_count=unmapped_reads,

            metadata=sample_metadata.get(sample_id, {}),
        )
        results.append(sample_result)

        # Print summary
        hdr_pct = summary.get('hdr_total_rate', 0) * 100
        nhej_pct = summary.get('nhej_mmej_total_rate', 0) * 100
        click.echo(f"  Reads: {aligned_reads}, HDR: {hdr_pct:.1f}%, NHEJ/MMEJ: {nhej_pct:.1f}%")

    # Write output
    output_file = output_path / "reclassified_editing_outcomes.tsv"
    write_results_tsv(results, output_file)

    # Generate summary report
    summary_file = output_path / "summary_report.md"
    generate_summary_report(results, summary_file)

    click.echo(f"\nReclassification complete!")
    click.echo(f"Processed {len(results)} samples")
    click.echo(f"Results written to: {output_file}")


@cli.command('aggregate')
@click.option('--input', '-i', 'input_file', type=click.Path(exists=True), required=True,
              help='Per-sample results TSV from TRACE pipeline')
@click.option('--output', '-o', type=click.Path(), required=True,
              help='Output unified comparison table TSV')
@click.option('--bio-sample-col', type=str, default='bio_sample',
              help='Column name for biological sample grouping (default: bio_sample)')
@click.option('--min-reads', type=int, default=1000,
              help='Minimum reads threshold for QC (default: 1000)')
@click.option('--include-crispresso/--no-crispresso', default=True,
              help='Include CRISPResso comparison columns (default: enabled)')
@click.option('--metadata-cols', type=str, default=None,
              help='Comma-separated list of additional metadata columns to include')
def aggregate(input_file, output, bio_sample_col, min_reads, include_crispresso, metadata_cols):
    """
    Generate unified comparison table with replicate aggregation.

    Takes per-sample results from TRACE pipeline, groups by biological sample,
    calculates mean/SEM for all metrics, and applies QC flags based on read counts.

    QC Flags:
      - good: All replicates have >= min_reads
      - all_low_reads: All replicates have < min_reads (still used)
      - one_low_read_removed: Low-read replicates excluded from mean
      - no_data: No data available

    \\b
    Example:
      trace aggregate \\
        --input results/per_sample_editing_outcomes.tsv \\
        --output results/unified_comparison_table.tsv \\
        --bio-sample-col bio_sample \\
        --min-reads 1000
    """
    import logging
    import pandas as pd

    from .analysis.aggregation import (
        generate_unified_comparison_table,
        write_unified_comparison_table,
    )

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

    # Load per-sample results
    click.echo(f"Loading per-sample results from: {input_file}")
    try:
        df = pd.read_csv(input_file, sep='\t')
    except Exception as e:
        click.echo(f"Error loading input file: {e}", err=True)
        sys.exit(1)

    click.echo(f"  Loaded {len(df)} samples")

    # Check for required columns
    if bio_sample_col not in df.columns:
        click.echo(f"Error: Column '{bio_sample_col}' not found in input file", err=True)
        click.echo(f"Available columns: {list(df.columns)[:10]}...", err=True)
        sys.exit(1)

    # Parse metadata columns
    metadata_list = None
    if metadata_cols:
        metadata_list = [c.strip() for c in metadata_cols.split(',')]

    # Generate unified comparison table
    click.echo(f"Generating unified comparison table...")
    click.echo(f"  Grouping by: {bio_sample_col}")
    click.echo(f"  Min reads threshold: {min_reads}")

    try:
        unified_df = generate_unified_comparison_table(
            df=df,
            bio_sample_col=bio_sample_col,
            min_reads=min_reads,
            include_crispresso=include_crispresso,
            metadata_cols=metadata_list,
        )
    except Exception as e:
        click.echo(f"Error generating unified table: {e}", err=True)
        sys.exit(1)

    # Write output
    output_path = Path(output)
    write_unified_comparison_table(unified_df, output_path)

    # Print summary
    n_bio_samples = len(unified_df)
    n_good = (unified_df['quality_flag'] == 'good').sum()
    n_low_reads = (unified_df['quality_flag'] == 'all_low_reads').sum()
    n_partial = (unified_df['quality_flag'] == 'one_low_read_removed').sum()

    click.echo(f"\n--- Summary ---")
    click.echo(f"Biological samples: {n_bio_samples}")
    click.echo(f"  Good (all replicates >= {min_reads} reads): {n_good}")
    click.echo(f"  All low reads: {n_low_reads}")
    click.echo(f"  Partial (some low replicates removed): {n_partial}")

    # Show mean HDR rate
    if 'dedup_hdr_total_pct_mean' in unified_df.columns:
        mean_hdr = unified_df['dedup_hdr_total_pct_mean'].mean()
        click.echo(f"\nMean HDR rate: {mean_hdr:.2f}%")

    click.echo(f"\nUnified comparison table written to: {output_path}")


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
