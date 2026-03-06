#!/usr/bin/env python3
"""
Classify reads from a single sample using TRACE.

Called by Snakemake with snakemake object providing:
- snakemake.input: Input files (BAM)
- snakemake.output: Output files (classification TSV)
- snakemake.params: Parameters (reference, guide, hdr_template, etc.)
- snakemake.log: Log file path
- snakemake.wildcards: Sample ID
"""

import sys
from pathlib import Path
import pandas as pd
import pysam

# Import TRACE modules
try:
    from trace_crispr.config import LocusConfig, NucleaseType, parse_sequence_input
    from trace_crispr.core.classification import classify_read, get_hdr_signature_positions
except ModuleNotFoundError as e:
    # Write error to log and exit
    log_path = Path(snakemake.log[0])
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text(f"ERROR: Failed to import trace_crispr: {e}\n")
    sys.exit(1)


def main():
    # Get parameters
    sample_id = snakemake.wildcards.sample
    reference = snakemake.params.reference
    guide = snakemake.params.guide
    hdr_template = snakemake.params.hdr_template

    analysis_window = snakemake.params.get('analysis_window', 10)
    large_deletion_min = snakemake.params.get('large_deletion_min', 50)
    hdr_threshold = snakemake.params.get('hdr_threshold', 0.8)

    log_path = Path(snakemake.log[0])
    log_path.parent.mkdir(parents=True, exist_ok=True)

    # Check if sample has guide/donor
    if not guide or not hdr_template or pd.isna(guide) or pd.isna(hdr_template):
        # Create empty classification file
        Path(snakemake.output.tsv).write_text("read_name\toutcome\thdr_fraction\n")
        log_path.write_text(f"Skipped {sample_id} - no guide or donor\n")
        return

    try:
        # Parse sequences
        ref_seq = parse_sequence_input(reference)
        donor_seq = parse_sequence_input(hdr_template)
        guide_seq = guide.upper()

        # Create locus config
        locus = LocusConfig(
            name=f"{sample_id}_locus",
            reference=ref_seq,
            hdr_template=donor_seq,
            guide=guide_seq,
            nuclease=NucleaseType.CAS9
        )
        locus.analyze()

        # Get HDR signature positions - must align template to reference first
        template_offset = locus.template_offset or 0
        aligned_ref = ref_seq[template_offset:template_offset + len(donor_seq)]
        hdr_signature = get_hdr_signature_positions(aligned_ref, donor_seq)

        # Get cut site from locus
        cut_site = locus.guide_info.cleavage_site if locus.guide_info else 0

        # Triple-aligner approach: Try aligners in order and use first successful
        # Order: BWA, BBMap, minimap2 (matching core TRACE behavior)
        aligner_paths = [
            ('bwa', snakemake.input.bwa),
            ('bbmap', snakemake.input.bbmap),
            ('minimap2', snakemake.input.minimap2)
        ]

        results = []
        aligner_used = None

        for aligner_name, bam_path in aligner_paths:
            bam_file = Path(bam_path)
            if not bam_file.exists():
                continue

            try:
                with pysam.AlignmentFile(str(bam_path), "rb") as bam:
                    # Test if BAM is readable and has alignments
                    for read in bam:
                        if read.is_unmapped:
                            continue

                        result = classify_read(
                            read,
                            hdr_signature,
                            cut_site,
                            window_size=analysis_window,
                            large_del_min_size=large_deletion_min,
                            hdr_threshold=hdr_threshold
                        )

                        results.append({
                            'read_name': read.query_name,
                            'outcome': result.outcome.name,
                            'hdr_match_fraction': result.details.get('hdr_match_fraction', 0.0)
                        })

                aligner_used = aligner_name
                break  # Successfully classified with this aligner

            except Exception as e:
                # This aligner failed, try next one
                log_path.write_text(f"Warning: {aligner_name} failed ({e}), trying next aligner...\n")
                continue

        if aligner_used is None:
            raise ValueError("All aligners failed - no successful alignments found")

        # Write results
        df = pd.DataFrame(results)
        df.to_csv(snakemake.output.tsv, sep='\t', index=False)

        log_path.write_text(f"Classified {len(results)} reads for {sample_id} using {aligner_used} aligner\n")

    except Exception as e:
        # Log error and create empty file
        import traceback
        log_path.write_text(f"ERROR classifying {sample_id}: {e}\n{traceback.format_exc()}\n")
        Path(snakemake.output.tsv).write_text("read_name\toutcome\thdr_fraction\n")
        raise


if __name__ == '__main__':
    main()
