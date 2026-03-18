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

# Add TRACE package to path (for script: directive execution)
TRACE_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(TRACE_ROOT))

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


def reverse_complement(seq):
    """Return reverse complement of a DNA sequence."""
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
            'a': 't', 't': 'a', 'g': 'c', 'c': 'g', 'N': 'N', 'n': 'n'}
    return ''.join(comp.get(b, b) for b in reversed(seq))


def detect_donor_capture(read, donor_seq, ref_seq, cut_site, min_insertion_size=10, min_match_fraction=0.8):
    """Detect if a read has an insertion containing donor sequence (NHEJ-mediated donor capture).

    This detects cases where donor DNA was captured at the DSB via NHEJ rather than
    proper HDR. The signature is a large insertion near the cut site that matches
    the donor template sequence.

    Args:
        read: pysam AlignedSegment
        donor_seq: Donor template sequence (correctly oriented)
        ref_seq: Reference sequence
        cut_site: Position of Cas9 cut site in reference coordinates
        min_insertion_size: Minimum insertion size to consider (default 10bp)
        min_match_fraction: Minimum fraction of insertion matching donor (default 0.8)

    Returns:
        dict with keys:
            - is_donor_capture: bool
            - insertion_size: int (size of largest insertion)
            - insertion_seq: str (sequence of insertion)
            - donor_match_start: int (where in donor the insertion matches)
            - donor_match_end: int
            - microhomology_5p: str (microhomology at 5' junction)
            - microhomology_3p: str (microhomology at 3' junction)
        or None if no donor capture detected
    """
    if read.is_unmapped or read.cigartuples is None:
        return None

    read_seq = read.query_sequence
    if not read_seq:
        return None

    donor_upper = donor_seq.upper()

    # Parse CIGAR to find insertions
    ref_pos = read.reference_start
    read_pos = 0

    insertions_near_cut = []

    for op, length in read.cigartuples:
        if op == 0:  # M (match/mismatch)
            ref_pos += length
            read_pos += length
        elif op == 1:  # I (insertion)
            # Check if insertion is near cut site
            distance_to_cut = abs(ref_pos - cut_site)
            if distance_to_cut <= 50 and length >= min_insertion_size:
                insertion_seq = read_seq[read_pos:read_pos + length]
                insertions_near_cut.append({
                    'size': length,
                    'seq': insertion_seq,
                    'ref_pos': ref_pos,
                    'read_pos': read_pos,
                    'distance_to_cut': distance_to_cut
                })
            read_pos += length
        elif op == 2:  # D (deletion)
            ref_pos += length
        elif op == 4:  # S (soft clip)
            read_pos += length
        elif op == 5:  # H (hard clip)
            pass

    if not insertions_near_cut:
        return None

    # Check if any insertion matches the donor
    for ins in insertions_near_cut:
        ins_seq = ins['seq'].upper()

        # Search for insertion sequence in donor
        # Allow for partial matches (the donor may have been trimmed)
        best_match_start = -1
        best_match_len = 0

        # Try to find the insertion sequence in the donor
        for start in range(len(donor_upper) - min_insertion_size + 1):
            # Try different match lengths
            for match_len in range(min(len(ins_seq), len(donor_upper) - start), min_insertion_size - 1, -1):
                donor_region = donor_upper[start:start + match_len]
                ins_region = ins_seq[:match_len]

                # Count matches
                matches = sum(1 for a, b in zip(donor_region, ins_region) if a == b)
                match_frac = matches / match_len

                if match_frac >= min_match_fraction and match_len > best_match_len:
                    best_match_start = start
                    best_match_len = match_len

        if best_match_start >= 0 and best_match_len >= min_insertion_size:
            # Found donor capture!
            # Look for microhomology at junctions

            # 5' junction: ref before insertion vs donor start
            ref_before = ref_seq[max(0, ins['ref_pos']-10):ins['ref_pos']].upper()
            donor_start_region = donor_upper[max(0, best_match_start-10):best_match_start+5].upper()
            microhomology_5p = find_microhomology(ref_before, donor_upper[best_match_start:best_match_start+10])

            # 3' junction: donor end vs ref after insertion
            donor_end = donor_upper[best_match_start + best_match_len - 10:best_match_start + best_match_len].upper()
            ref_after = ref_seq[ins['ref_pos']:ins['ref_pos']+15].upper()
            microhomology_3p = find_microhomology(donor_end, ref_after)

            return {
                'is_donor_capture': True,
                'insertion_size': ins['size'],
                'insertion_seq': ins['seq'],
                'donor_match_start': best_match_start,
                'donor_match_end': best_match_start + best_match_len,
                'match_length': best_match_len,
                'distance_to_cut': ins['distance_to_cut'],
                'microhomology_5p': microhomology_5p,
                'microhomology_3p': microhomology_3p
            }

    return None


def find_microhomology(seq1, seq2, min_length=2, max_length=10):
    """Find microhomology between end of seq1 and start of seq2.

    Returns the longest microhomology sequence found, or empty string if none.
    """
    best_homology = ""

    for length in range(min(max_length, len(seq1), len(seq2)), min_length - 1, -1):
        if seq1[-length:] == seq2[:length]:
            return seq1[-length:]

    return best_homology


def find_best_alignment(ref_seq, donor_seq):
    """Find best alignment offset and orientation for donor in reference.

    Returns:
        (offset, donor_to_use, is_rc): Best offset, donor sequence to use, and whether it's RC
    """
    ref_upper = ref_seq.upper()
    donor_fwd = donor_seq
    donor_rc = reverse_complement(donor_seq)

    best_offset = 0
    best_matches = 0
    best_is_rc = False

    # Try forward orientation
    for offset in range(max(1, len(ref_upper) - len(donor_fwd) + 1)):
        aligned = ref_upper[offset:offset + len(donor_fwd)]
        if len(aligned) == len(donor_fwd):
            matches = sum(1 for a, b in zip(aligned, donor_fwd.upper()) if a == b)
            if matches > best_matches:
                best_matches = matches
                best_offset = offset
                best_is_rc = False

    # Try reverse complement
    for offset in range(max(1, len(ref_upper) - len(donor_rc) + 1)):
        aligned = ref_upper[offset:offset + len(donor_rc)]
        if len(aligned) == len(donor_rc):
            matches = sum(1 for a, b in zip(aligned, donor_rc.upper()) if a == b)
            if matches > best_matches:
                best_matches = matches
                best_offset = offset
                best_is_rc = True

    return best_offset, donor_rc if best_is_rc else donor_fwd, best_is_rc


def find_edit_region(ref_aligned, donor_aligned, cut_site=None, max_gap=30):
    """Find the edit region by clustering SNVs near the cut site.

    Intended edits are clustered near the cut site. Flanking sequences
    (SSO14 motifs, HUHe sites) create SNVs far from the cut site.

    Args:
        ref_aligned: Reference sequence at alignment position
        donor_aligned: Donor sequence (already oriented correctly)
        cut_site: Position of the cut site (if known)
        max_gap: Maximum gap between SNVs to be in same cluster

    Returns:
        (edit_start, edit_end): Boundaries of the edit region
        - Only SNVs within this region are counted as intended edits
    """
    if len(ref_aligned) != len(donor_aligned):
        return 0, len(donor_aligned)

    n = len(donor_aligned)

    # Find all SNV positions
    snv_positions = []
    for i in range(n):
        if ref_aligned[i].upper() != donor_aligned[i].upper():
            snv_positions.append(i)

    if not snv_positions:
        return 0, n

    # If only one cluster of SNVs, use the entire range
    if len(snv_positions) == 1:
        return snv_positions[0], snv_positions[0] + 1

    # Find clusters by looking for gaps > max_gap between consecutive SNVs
    clusters = []
    cluster_start = snv_positions[0]
    cluster_end = snv_positions[0]

    for i in range(1, len(snv_positions)):
        if snv_positions[i] - snv_positions[i-1] > max_gap:
            # Gap found - save current cluster and start new one
            clusters.append((cluster_start, cluster_end + 1))
            cluster_start = snv_positions[i]
        cluster_end = snv_positions[i]

    # Don't forget the last cluster
    clusters.append((cluster_start, cluster_end + 1))

    if len(clusters) == 1:
        # Only one cluster - all SNVs are intended edits
        return clusters[0]

    # Multiple clusters - choose the one nearest the center of the donor
    # (edit site is typically in the middle, flanking sequences at ends)
    center = n // 2

    if cut_site is not None:
        # If cut site is known, use that as the reference point
        center = cut_site

    # Find cluster closest to center
    best_cluster = clusters[0]
    best_distance = abs((clusters[0][0] + clusters[0][1]) / 2 - center)

    for cluster in clusters[1:]:
        cluster_center = (cluster[0] + cluster[1]) / 2
        distance = abs(cluster_center - center)
        if distance < best_distance:
            best_distance = distance
            best_cluster = cluster

    return best_cluster


def main():
    # Get parameters
    sample_id = snakemake.wildcards.sample
    reference = snakemake.params.reference
    guide = snakemake.params.guide
    hdr_template = snakemake.params.hdr_template
    # Contamination check template: used when sample lacks HDR template but we want to check for HDR contamination
    contamination_check_template = snakemake.params.get('contamination_check_template', '')

    analysis_window = snakemake.params.get('analysis_window', 10)
    large_deletion_min = snakemake.params.get('large_deletion_min', 50)
    hdr_threshold = snakemake.params.get('hdr_threshold', 0.8)

    log_path = Path(snakemake.log[0])
    log_path.parent.mkdir(parents=True, exist_ok=True)

    # Check if sample has guide - guide is required for classification
    # HDR template is optional: without it, we still detect NHEJ vs WT
    if not guide or pd.isna(guide):
        # No guide = no DSB expected, skip classification
        Path(snakemake.output.tsv).write_text("read_name\toutcome\thdr_fraction\n")
        Path(snakemake.output.hdr_detail).write_text("read_name\tsnv_position\tsnv_detected\n")
        log_path.write_text(f"Skipped {sample_id} - no guide sequence\n")
        return

    # Check if we have HDR template for HDR detection
    has_hdr_template = hdr_template and not pd.isna(hdr_template) and str(hdr_template).strip()

    # Contamination check mode: no sample-specific HDR template, but use default to check for contamination/swap
    contamination_check_mode = False
    effective_hdr_template = hdr_template  # The template actually used for analysis
    if not has_hdr_template and contamination_check_template and str(contamination_check_template).strip():
        contamination_check_mode = True
        effective_hdr_template = contamination_check_template
        has_hdr_template = True  # Use contamination template for HDR detection

    try:
        # Parse sequences (converts to uppercase)
        ref_seq = parse_sequence_input(reference)
        guide_seq = guide.upper()

        # Parse HDR template if available
        donor_seq = None
        if has_hdr_template:
            donor_seq = parse_sequence_input(effective_hdr_template)

        # Variables for HDR detection (only used if we have HDR template)
        template_offset = 0
        donor_aligned = None
        is_rc = False
        hdr_signature = []
        flanking_signatures = []
        aligned_ref = ref_seq
        edit_start, edit_end = 0, 0

        if has_hdr_template:
            # Find best alignment - handles both forward and reverse complement donors
            template_offset, donor_aligned, is_rc = find_best_alignment(ref_seq, donor_seq)

            # Create locus config with correctly oriented donor
            locus = LocusConfig(
                name=f"{sample_id}_locus",
                reference=ref_seq,
                hdr_template=donor_aligned,
                guide=guide_seq,
                nuclease=NucleaseType.CAS9
            )
            locus.analyze()

            # Get aligned reference region for comparison
            aligned_ref = ref_seq[template_offset:template_offset + len(donor_aligned)]

            # Get cut site from locus (relative to aligned region)
            cut_site = locus.guide_info.cleavage_site if locus.guide_info else len(donor_aligned) // 2
            cut_site_in_donor = cut_site - template_offset  # Convert to donor-relative position

            # Get all signature positions via de novo comparison
            all_signatures = get_hdr_signature_positions(aligned_ref, donor_aligned)

            # Find the edit region by clustering SNVs near the cut site
            # Intended edits cluster near cut site; flanking sequences (SSO14/HUHe) are at donor ends
            edit_start, edit_end = find_edit_region(aligned_ref, donor_aligned, cut_site_in_donor)

            # Filter to only include signatures in the edit region
            hdr_signature = [(pos, wt, hdr) for pos, wt, hdr in all_signatures
                             if edit_start <= pos < edit_end]

            # Flanking signatures are SNVs outside the edit region
            # These can be checked to detect unintended incorporation of flanking sequences
            flanking_signatures = [(pos, wt, hdr) for pos, wt, hdr in all_signatures
                                   if pos < edit_start or pos >= edit_end]
        else:
            # No HDR template - find cut site from guide position in reference
            guide_pos = ref_seq.upper().find(guide_seq.upper())
            if guide_pos == -1:
                # Try reverse complement
                guide_rc = reverse_complement(guide_seq).upper()
                guide_pos = ref_seq.upper().find(guide_rc)

            if guide_pos >= 0:
                # Cut site is 3bp upstream of PAM (end of guide) for Cas9
                cut_site = guide_pos + len(guide_seq) - 3
            else:
                # Guide not found, use middle of reference
                cut_site = len(ref_seq) // 2

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
                            ref_offset=template_offset,  # Critical: convert signature positions to reference coords
                            window_size=analysis_window,
                            large_del_min_size=large_deletion_min,
                            hdr_threshold=hdr_threshold
                        )

                        # Check for unintended incorporation of flanking sequences
                        flanking_matches = 0
                        flanking_covered = 0
                        if flanking_signatures:
                            from trace_crispr.core.classification import count_hdr_signature_matches
                            flanking_matches, flanking_covered, _ = count_hdr_signature_matches(
                                read, flanking_signatures, ref_offset=template_offset
                            )

                        # Check for donor capture (NHEJ-mediated insertion of donor fragment)
                        # Only check if we have HDR template
                        donor_capture = None
                        if has_hdr_template and donor_aligned:
                            donor_capture = detect_donor_capture(
                                read, donor_aligned, ref_seq, cut_site,
                                min_insertion_size=15, min_match_fraction=0.85
                            )

                        # Determine final outcome
                        outcome = result.outcome.name
                        is_donor_capture = False
                        donor_capture_size = 0

                        if donor_capture and donor_capture['is_donor_capture']:
                            # Override outcome to DONOR_CAPTURE
                            is_donor_capture = True
                            donor_capture_size = donor_capture['insertion_size']
                            outcome = 'DONOR_CAPTURE'

                        # In contamination check mode, HDR detection means possible contamination/sample swap
                        if contamination_check_mode and outcome in ('HDR_PERFECT', 'HDR_IMPERFECT'):
                            outcome = 'HDR_CONTAMINATION'  # Flag as potential contamination

                        results.append({
                            'read_name': read.query_name,
                            'outcome': outcome,
                            'hdr_match_fraction': result.details.get('hdr_match_fraction', 0.0),
                            'flanking_matches': flanking_matches,
                            'flanking_covered': flanking_covered,
                            'is_donor_capture': is_donor_capture,
                            'donor_capture_size': donor_capture_size
                        })

                aligner_used = aligner_name
                break  # Successfully classified with this aligner

            except Exception as e:
                # This aligner failed, try next one
                log_path.write_text(f"Warning: {aligner_name} failed ({e}), trying next aligner...\n")
                continue

        if aligner_used is None:
            raise ValueError("All aligners failed - no successful alignments found")

        # Deduplicate paired-end reads: keep the read with highest hdr_match_fraction per read pair
        # This handles cases where R1 covers HDR region but R2 doesn't (or vice versa)
        df = pd.DataFrame(results)
        original_count = len(df)

        if len(df) > 0:
            # Sort by hdr_match_fraction descending, then keep first occurrence of each read name
            df = df.sort_values('hdr_match_fraction', ascending=False)
            df = df.drop_duplicates(subset='read_name', keep='first')
            df = df.sort_index()  # Restore original order

        dedup_count = len(df)
        df.to_csv(snakemake.output.tsv, sep='\t', index=False)

        # Write HDR SNV detail file (placeholder for now)
        Path(snakemake.output.hdr_detail).write_text("read_name\tsnv_position\tsnv_detected\n")

        # Build log notes based on whether HDR template was available
        if contamination_check_mode:
            rc_note = " [donor RC]" if is_rc else ""
            edit_note = " [CONTAMINATION CHECK MODE - checking for HDR from other samples]"
            excluded_note = f" ({len(all_signatures) - len(hdr_signature)} flanking SNVs excluded)" if len(all_signatures) > len(hdr_signature) else ""
        elif has_hdr_template:
            rc_note = " [donor RC]" if is_rc else ""
            edit_note = f" [edit region: {edit_start}-{edit_end} of {len(donor_aligned)}bp]"
            excluded_note = f" ({len(all_signatures) - len(hdr_signature)} flanking SNVs excluded)" if len(all_signatures) > len(hdr_signature) else ""
        else:
            rc_note = ""
            edit_note = " [NO HDR TEMPLATE - NHEJ/WT detection only]"
            excluded_note = ""
            all_signatures = []

        # Calculate flanking incorporation stats
        flanking_note = ""
        if len(df) > 0 and 'flanking_matches' in df.columns and flanking_signatures:
            reads_with_flanking = (df['flanking_matches'] > 0).sum()
            if reads_with_flanking > 0:
                flanking_pct = 100 * reads_with_flanking / len(df)
                flanking_note = f"\n  WARNING: {reads_with_flanking} reads ({flanking_pct:.2f}%) show unintended flanking sequence incorporation"

        # Calculate donor capture stats
        donor_capture_note = ""
        if len(df) > 0 and 'is_donor_capture' in df.columns:
            donor_captures = df['is_donor_capture'].sum()
            if donor_captures > 0:
                capture_pct = 100 * donor_captures / len(df)
                donor_capture_note = f"\n  NOTE: {donor_captures} reads ({capture_pct:.2f}%) show NHEJ-mediated donor capture (not true HDR)"

        # Check for contamination (HDR in no-donor sample)
        contamination_note = ""
        if contamination_check_mode and len(df) > 0:
            contamination_reads = (df['outcome'] == 'HDR_CONTAMINATION').sum()
            if contamination_reads > 0:
                contamination_pct = 100 * contamination_reads / len(df)
                contamination_note = f"\n  ⚠️  WARNING: {contamination_reads} reads ({contamination_pct:.2f}%) show HDR edits - possible sample contamination or swap!"

        log_path.write_text(f"Classified {dedup_count} read pairs for {sample_id} using {aligner_used} aligner{rc_note}{edit_note} ({len(hdr_signature)} signature positions{excluded_note}, deduplicated from {original_count} alignments){flanking_note}{donor_capture_note}{contamination_note}\n")

    except Exception as e:
        # Log error and create empty files
        import traceback
        log_path.write_text(f"ERROR classifying {sample_id}: {e}\n{traceback.format_exc()}\n")
        Path(snakemake.output.tsv).write_text("read_name\toutcome\thdr_fraction\n")
        Path(snakemake.output.hdr_detail).write_text("read_name\tsnv_position\tsnv_detected\n")
        raise


if __name__ == '__main__':
    main()
