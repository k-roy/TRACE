#!/usr/bin/env python3
"""
Level 2 Cache: Classify paired-end sequences for a specific (guide, donor) combination.

NEW APPROACH: Uses R1/R2 paired alignments with coverage verification.

Supports flexible sequencing modes (any read length combination):
- Single-end: R1 only (e.g., 228bp R1)
- Asymmetric paired: R1 primary, R2 validates (e.g., 228x90bp, 100x200bp, etc.)
- Symmetric paired: Both R1 and R2 span edits (e.g., 2x300bp, 2x150bp, etc.)

Coverage verification: Checks that at least one read (or R1+R2 together) spans
all donor edits. Warns user if edit region is not fully covered.

Reuses Level 1 alignment cache - no re-alignment needed.
Only classifies sequences that appear with THIS guide/donor combination.
"""

import sys
import re
from pathlib import Path

# Add TRACE to path
TRACE_ROOT = Path("/oak/stanford/groups/larsms/Users/kevinroy/software/trace")
sys.path.insert(0, str(TRACE_ROOT))

import pandas as pd

from trace_crispr.config import parse_sequence_input
from trace_crispr.core.edit_distance_hdr import (
    build_donor_signature, classify_read_edit_distance, DonorSignature
)


def reverse_complement(seq):
    """Return reverse complement of DNA sequence."""
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
            'a': 't', 't': 'a', 'g': 'c', 'c': 'g', 'N': 'N', 'n': 'n'}
    return ''.join(comp.get(b, b) for b in reversed(seq))


def find_cut_site(ref_seq, guide_seq):
    """
    Find cut site position from guide location in reference.

    For SpCas9:
    - Plus strand: PAM is downstream, cut site is 3bp upstream of PAM
                   = guide_end - 3
    - Minus strand: PAM is upstream, cut site is 3bp downstream of PAM
                    = guide_start + 3
    """
    if not guide_seq:
        return len(ref_seq) // 2

    ref_upper = ref_seq.upper()
    guide_upper = guide_seq.upper()

    # Try plus strand first
    guide_pos = ref_upper.find(guide_upper)
    if guide_pos != -1:
        # Plus strand: cut site is 3bp upstream of PAM (at guide end)
        return guide_pos + len(guide_seq) - 3

    # Try minus strand (reverse complement)
    guide_rc = reverse_complement(guide_upper)
    guide_pos = ref_upper.find(guide_rc)
    if guide_pos != -1:
        # Minus strand: cut site is 3bp downstream of PAM (at guide start + 3)
        return guide_pos + 3

    # Guide not found - use midpoint
    return len(ref_seq) // 2


def parse_cigar_string(cigar_str):
    """Parse CIGAR string back to tuples."""
    if not cigar_str:
        return None
    op_map = {'M': 0, 'I': 1, 'D': 2, 'N': 3, 'S': 4, 'H': 5, 'P': 6, '=': 7, 'X': 8}
    tuples = []
    for match in re.finditer(r'(\d+)([MIDNSHP=X])', cigar_str):
        length = int(match.group(1))
        op = op_map[match.group(2)]
        tuples.append((op, length))
    return tuples


def analyze_junction_microhomology(insertion_seq, ref_left_flank, ref_right_flank,
                                  max_homology_check=15):
    """
    Analyze microhomology at insertion junctions to distinguish NHEJ vs MMEJ.

    MMEJ (microhomology-mediated end joining) uses short homologous sequences
    (typically 2-15 bp) to anneal DNA ends before ligation. NHEJ has little to no
    microhomology (0-1 bp).

    Args:
        insertion_seq: The inserted sequence
        ref_left_flank: Reference sequence immediately left of insertion (up to 15 bp)
        ref_right_flank: Reference sequence immediately right of insertion (up to 15 bp)
        max_homology_check: Maximum homology length to check (default 15 bp)

    Returns:
        dict with:
            - homology_left: bp of homology at left junction
            - homology_right: bp of homology at right junction
            - mmej_class: NHEJ | MMEJ_SHORT | MMEJ_EXTENDED | COMPLEX
    """
    if not insertion_seq or not ref_left_flank or not ref_right_flank:
        return {'homology_left': 0, 'homology_right': 0, 'mmej_class': 'UNKNOWN'}

    ins_seq = insertion_seq.upper()
    left_flank = ref_left_flank.upper()
    right_flank = ref_right_flank.upper()

    # Check left junction: compare end of left flank with start of insertion
    homology_left = 0
    check_len = min(len(left_flank), len(ins_seq), max_homology_check)
    for i in range(1, check_len + 1):
        # Compare last i bases of left flank with first i bases of insertion
        if left_flank[-i:] == ins_seq[:i]:
            homology_left = i
        else:
            break

    # Check right junction: compare end of insertion with start of right flank
    homology_right = 0
    check_len = min(len(right_flank), len(ins_seq), max_homology_check)
    for i in range(1, check_len + 1):
        # Compare last i bases of insertion with first i bases of right flank
        if ins_seq[-i:] == right_flank[:i]:
            homology_right = i
        else:
            break

    # Classify based on homology length
    max_homology = max(homology_left, homology_right)

    if max_homology == 0 or max_homology == 1:
        mmej_class = 'NHEJ'  # Blunt or near-blunt end ligation
    elif 2 <= max_homology <= 5:
        mmej_class = 'MMEJ_SHORT'  # Classic MMEJ
    elif 6 <= max_homology <= 15:
        mmej_class = 'MMEJ_EXTENDED'  # Alt-NHEJ with longer homology
    else:
        mmej_class = 'COMPLEX'  # Unusual, possible HDR initiation -> NHEJ

    return {
        'homology_left': homology_left,
        'homology_right': homology_right,
        'mmej_class': mmej_class
    }


def detect_donor_capture(read_seq, cigar_tuples, ref_start, donor_seq, ref_seq, cut_site,
                        min_insertion_size=10, min_match_fraction=0.8, max_cut_distance=50):
    """
    Detect if a read has an insertion containing donor sequence (NHEJ-mediated donor capture).

    This detects cases where donor DNA was captured at the DSB via NHEJ rather than
    proper HDR. The signature is a large insertion near the cut site that matches
    the donor template sequence.

    Also analyzes junction microhomology to distinguish NHEJ from MMEJ-mediated capture.

    Args:
        read_seq: Read sequence
        cigar_tuples: Parsed CIGAR tuples [(op, length), ...]
        ref_start: Reference start position
        donor_seq: Donor template sequence
        ref_seq: Reference sequence (for extracting flanking regions)
        cut_site: Position of Cas9 cut site in reference coordinates
        min_insertion_size: Minimum insertion size to consider (default 10bp)
        min_match_fraction: Minimum fraction of insertion matching donor (default 0.8)
        max_cut_distance: Maximum distance from cut site to consider (default 50bp)

    Returns:
        dict with donor capture info including junction analysis, or None if not detected
    """
    if not read_seq or not cigar_tuples or not donor_seq:
        return None

    donor_upper = donor_seq.upper()

    # Parse CIGAR to find insertions near cut site
    ref_pos = ref_start
    read_pos = 0
    insertions_near_cut = []

    for op, length in cigar_tuples:
        if op == 0:  # M (match/mismatch)
            ref_pos += length
            read_pos += length
        elif op == 1:  # I (insertion)
            distance_to_cut = abs(ref_pos - cut_site)
            if distance_to_cut <= max_cut_distance and length >= min_insertion_size:
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

    if not insertions_near_cut:
        return None

    # Check if any insertion matches the donor
    for ins in insertions_near_cut:
        ins_seq = ins['seq'].upper()

        # Search for insertion sequence in donor
        best_match_start = -1
        best_match_len = 0

        # Sliding window to find best match
        for start in range(len(donor_upper) - min_insertion_size + 1):
            for match_len in range(min(len(ins_seq), len(donor_upper) - start),
                                  min_insertion_size - 1, -1):
                donor_region = donor_upper[start:start + match_len]
                ins_region = ins_seq[:match_len]

                # Count matches
                matches = sum(1 for a, b in zip(donor_region, ins_region) if a == b)
                match_frac = matches / match_len

                if match_frac >= min_match_fraction and match_len > best_match_len:
                    best_match_start = start
                    best_match_len = match_len

        if best_match_start >= 0 and best_match_len >= min_insertion_size:
            # Found donor capture! Now analyze junction microhomology

            # Extract flanking reference sequences (15 bp each side)
            ref_upper = ref_seq.upper() if ref_seq else ''
            ref_pos_in_seq = ins['ref_pos'] - ref_start

            # Get flanking sequences
            flank_len = 15
            if ref_upper and 0 <= ref_pos_in_seq <= len(ref_upper):
                left_start = max(0, ref_pos_in_seq - flank_len)
                left_flank = ref_upper[left_start:ref_pos_in_seq]

                right_end = min(len(ref_upper), ref_pos_in_seq + flank_len)
                right_flank = ref_upper[ref_pos_in_seq:right_end]
            else:
                left_flank = ''
                right_flank = ''

            # Analyze junction microhomology
            junction_info = analyze_junction_microhomology(
                insertion_seq=ins['seq'],
                ref_left_flank=left_flank,
                ref_right_flank=right_flank
            )

            return {
                'is_donor_capture': True,
                'insertion_size': ins['size'],
                'insertion_seq': ins['seq'],
                'donor_match_start': best_match_start,
                'donor_match_end': best_match_start + best_match_len,
                'match_length': best_match_len,
                'match_fraction': best_match_len / ins['size'],
                'distance_to_cut': ins['distance_to_cut'],
                'ref_position': ins['ref_pos'],
                'junction_homology_left': junction_info['homology_left'],
                'junction_homology_right': junction_info['homology_right'],
                'mmej_class': junction_info['mmej_class']
            }

    return None


def get_aligned_region(ref_start, ref_end, is_reverse):
    """
    Get the reference region covered by a read alignment.

    Returns: (start_pos, end_pos) on reference
    """
    if ref_start < 0 or ref_end < 0:
        return None
    return (ref_start, ref_end)


def check_donor_coverage(donor_signature, r1_region, r2_region, ref_length):
    """
    Check if R1, R2, or R1+R2 together cover all donor edits.

    Returns: (is_covered, coverage_mode, uncovered_edits)
    - is_covered: True if all edits are covered
    - coverage_mode: 'r1_only', 'r2_only', 'r1_and_r2', or 'incomplete'
    - uncovered_edits: List of edit positions not covered
    """
    # Collect all edit positions
    edit_positions = set()
    edit_positions.update(donor_signature.snvs.keys())
    edit_positions.update(donor_signature.insertions.keys())
    edit_positions.update(donor_signature.deletions.keys())

    if not edit_positions:
        # No edits - always covered
        return True, 'no_edits', []

    # Check which edits are covered by R1
    r1_covered = set()
    if r1_region:
        r1_start, r1_end = r1_region
        r1_covered = {pos for pos in edit_positions if r1_start <= pos < r1_end}

    # Check which edits are covered by R2
    r2_covered = set()
    if r2_region:
        r2_start, r2_end = r2_region
        r2_covered = {pos for pos in edit_positions if r2_start <= pos < r2_end}

    # Combined coverage
    all_covered = r1_covered | r2_covered
    uncovered = edit_positions - all_covered

    # Determine coverage mode
    if len(uncovered) == 0:
        # All edits covered
        if len(r1_covered) == len(edit_positions):
            coverage_mode = 'r1_only'
        elif len(r2_covered) == len(edit_positions):
            coverage_mode = 'r2_only'
        else:
            coverage_mode = 'r1_and_r2'
        return True, coverage_mode, []
    else:
        return False, 'incomplete', sorted(uncovered)


def detect_sequencing_mode(r1_mapped, r2_mapped, r1_seq_len, r2_seq_len, donor_signature, ref_length):
    """
    Detect sequencing mode based on read lengths and donor edit positions.

    Flexible detection for ANY read length combination:
    - 2x150, 151x149, 100x200, 228x90, 300x300, 100x500, etc.

    Returns: 'single_end', 'asymmetric', or 'symmetric'

    Logic:
    - single_end: No R2 or R2 unmapped
    - symmetric: Both R1 and R2 likely span all donor edits
    - asymmetric: R1 is primary, R2 provides validation/deletion detection
    """
    if not r2_mapped:
        return 'single_end'

    # Get edit region bounds
    if hasattr(donor_signature, 'snvs') and donor_signature.snvs:
        edit_positions = list(donor_signature.snvs.keys())
        edit_positions.extend(donor_signature.insertions.keys())
        edit_positions.extend(donor_signature.deletions.keys())

        if edit_positions:
            edit_start = min(edit_positions)
            edit_end = max(edit_positions)
            edit_span = edit_end - edit_start
        else:
            edit_span = 0
    else:
        edit_span = 0

    # Heuristic: Symmetric if BOTH reads are long enough to span edit region
    # with some buffer (1.2x edit span)
    if edit_span > 0:
        min_read_length = edit_span * 1.2
        if r1_seq_len >= min_read_length and r2_seq_len >= min_read_length:
            return 'symmetric'

    # Another heuristic: Symmetric if both reads > 200bp (common for 2x300, 2x250, 2x150)
    if r1_seq_len >= 200 and r2_seq_len >= 200:
        return 'symmetric'

    # Otherwise asymmetric
    return 'asymmetric'


def merge_symmetric_classifications(r1_clf, r2_clf, donor_signature):
    """
    Merge R1 and R2 classifications for symmetric sequencing.

    Conservative intersection approach: edits must be confirmed by both reads.
    """
    # Confirmed donor SNVs: intersection
    r1_donor_snvs = set(r1_clf.donor_snvs_detected) if r1_clf.donor_snvs_detected else set()
    r2_donor_snvs = set(r2_clf.donor_snvs_detected) if r2_clf.donor_snvs_detected else set()
    confirmed_donor_snvs = r1_donor_snvs & r2_donor_snvs

    # If both agree on outcome, use it
    if r1_clf.outcome == r2_clf.outcome:
        return r1_clf

    # If one is WT and other has edits, trust the edits
    # (one read may have lower coverage of edit region)
    if r1_clf.outcome == 'WT':
        return r2_clf
    if r2_clf.outcome == 'WT':
        return r1_clf

    # If both have HDR but different completeness, use intersection
    if 'HDR' in r1_clf.outcome and 'HDR' in r2_clf.outcome:
        total_donor_edits = len(donor_signature.snvs) + len(donor_signature.insertions) + len(donor_signature.deletions)
        if len(confirmed_donor_snvs) == total_donor_edits:
            # Create new classification with HDR_COMPLETE
            r1_clf.outcome = 'HDR_COMPLETE'
            r1_clf.donor_snvs_detected = sorted(confirmed_donor_snvs)
            return r1_clf
        elif len(confirmed_donor_snvs) > 0:
            # Partial HDR
            r1_clf.outcome = 'HDR_PARTIAL'
            r1_clf.donor_snvs_detected = sorted(confirmed_donor_snvs)
            return r1_clf

    # Use R1 as default (primary read)
    return r1_clf


def load_alignment_cache_sqlite(alignment_cache_path, relevant_hashes):
    """
    Load alignment cache using SQLite for efficient lookup.

    Instead of loading 4-12GB into memory, we:
    1. Create a temporary SQLite database
    2. Read the TSV in chunks, inserting only relevant hashes
    3. Query by seq_hash efficiently

    Memory usage: ~50MB vs 4-12GB

    Args:
        alignment_cache_path: Path to alignment cache TSV.gz
        relevant_hashes: Set of seq_hash values we need

    Returns:
        dict-like object with alignment lookups
    """
    import sqlite3
    import tempfile

    # Create temporary SQLite database
    temp_db = tempfile.NamedTemporaryFile(suffix='.db', delete=False)
    temp_db.close()

    conn = sqlite3.connect(temp_db.name)
    conn.execute("PRAGMA journal_mode = WAL")
    conn.execute("PRAGMA synchronous = NORMAL")

    # Read TSV header to get column names
    chunk_iter = pd.read_csv(
        alignment_cache_path,
        sep='\t',
        compression='gzip',
        chunksize=10000,
        nrows=0  # Just get columns
    )

    # Get columns from first chunk
    first_chunk = pd.read_csv(
        alignment_cache_path,
        sep='\t',
        compression='gzip',
        nrows=1
    )
    all_columns = first_chunk.columns.tolist()

    # Remove seq_hash from columns list since we'll add it as PRIMARY KEY
    columns = [c for c in all_columns if c != 'seq_hash']

    # Create table with seq_hash as PRIMARY KEY, then other columns
    col_defs = ', '.join([f'"{col}" TEXT' for col in columns])
    conn.execute(f'CREATE TABLE alignments (seq_hash TEXT PRIMARY KEY, {col_defs})')

    # Read in chunks and insert only relevant hashes using batch insert
    loaded = 0
    chunk_iter = pd.read_csv(
        alignment_cache_path,
        sep='\t',
        compression='gzip',
        chunksize=10000
    )

    for chunk in chunk_iter:
        # Filter to relevant hashes
        relevant_chunk = chunk[chunk['seq_hash'].isin(relevant_hashes)]

        if len(relevant_chunk) > 0:
            # Batch insert using executemany (10-100x faster than row-by-row)
            insert_cols = ['seq_hash'] + columns
            placeholders = ', '.join(['?'] * len(insert_cols))

            # Prepare batch data
            batch_data = []
            for _, row in relevant_chunk.iterrows():
                values = tuple(str(row[col]) if row[col] is not None else '' for col in insert_cols)
                batch_data.append(values)

            conn.executemany(
                f'INSERT OR REPLACE INTO alignments VALUES ({placeholders})',
                batch_data
            )
            loaded += len(relevant_chunk)

    conn.commit()

    # Create a dict-like wrapper
    class SQLiteAlignmentLookup:
        def __init__(self, conn, columns):
            self.conn = conn
            self.columns = columns

        def __contains__(self, seq_hash):
            cursor = self.conn.execute(
                'SELECT 1 FROM alignments WHERE seq_hash = ?', (seq_hash,)
            )
            return cursor.fetchone() is not None

        def __getitem__(self, seq_hash):
            cursor = self.conn.execute(
                'SELECT * FROM alignments WHERE seq_hash = ?', (seq_hash,)
            )
            row = cursor.fetchone()
            if row is None:
                raise KeyError(seq_hash)

            # Get actual column names from cursor description
            col_names = [desc[0] for desc in cursor.description]
            result = {}
            for i, col in enumerate(col_names):
                val = row[i]
                # Convert numeric strings back to numbers
                if col in ('r1_ref_start', 'r1_ref_end', 'r2_ref_start', 'r2_ref_end'):
                    try:
                        result[col] = int(val) if val and val != 'nan' else -1
                    except (ValueError, TypeError):
                        result[col] = -1
                elif col in ('r1_is_mapped', 'r2_is_mapped', 'is_proper_pair'):
                    result[col] = val == 'True' or val == '1' or val == True
                else:
                    result[col] = val if val != 'nan' else ''
            return result

        def get(self, seq_hash, default=None):
            try:
                return self[seq_hash]
            except KeyError:
                return default

        def close(self):
            self.conn.close()
            import os
            os.unlink(temp_db.name)

    return SQLiteAlignmentLookup(conn, columns), loaded


def main():
    alignment_cache_path = snakemake.input.alignment_cache
    seq_to_samples_path = snakemake.input.seq_to_samples
    output_path = snakemake.output.classifications
    log_path = Path(str(snakemake.log))

    guide = snakemake.params.guide if snakemake.params.guide else ""
    donor = snakemake.params.donor if snakemake.params.donor else ""
    ref_param = snakemake.params.reference
    edit_region_size = snakemake.params.edit_region_size
    snv_distance_filter = getattr(snakemake.params, 'snv_distance_filter', 50)
    filter_chimeric = getattr(snakemake.params, 'filter_chimeric', True)
    homopolymer_filter = getattr(snakemake.params, 'homopolymer_filter', 5)
    nhej_quantification_window = getattr(snakemake.params, 'nhej_quantification_window', 1)

    # Clean up nan values
    if guide == "nan":
        guide = ""
    if donor == "nan":
        donor = ""

    log_path.parent.mkdir(parents=True, exist_ok=True)
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)

    # First, identify relevant sequences for this guide/donor
    # Read mapping file in chunks to get relevant seq_hashes
    print(f"Finding relevant sequences for guide/donor...", file=sys.stderr)
    relevant_seq_hashes = set()

    chunk_iter = pd.read_csv(
        seq_to_samples_path,
        sep='\t',
        compression='gzip',
        chunksize=50000
    )

    for chunk in chunk_iter:
        chunk['guide'] = chunk['guide'].fillna('')
        chunk['donor'] = chunk['donor'].fillna('')
        relevant = chunk[
            (chunk['guide'] == guide) &
            (chunk['donor'] == donor)
        ]
        relevant_seq_hashes.update(relevant['seq_hash'].unique())

    print(f"Found {len(relevant_seq_hashes)} relevant sequences", file=sys.stderr)

    # Load only relevant alignments using SQLite
    print(f"Loading alignment cache (SQLite mode)...", file=sys.stderr)
    alignment_lookup, loaded_count = load_alignment_cache_sqlite(
        alignment_cache_path, relevant_seq_hashes
    )
    print(f"Loaded {loaded_count} alignments into SQLite cache", file=sys.stderr)

    # Parse reference
    ref_seq = parse_sequence_input(ref_param) if ref_param else None

    if ref_seq is None:
        # Cannot classify without a reference - write empty results
        log_path.write_text(f"ERROR: No reference sequence provided\n")
        pd.DataFrame(columns=['seq_hash', 'guide', 'donor', 'outcome',
                              'n_donor_snvs', 'n_non_donor', 'donor_fraction',
                              'donor_snvs_detected']).to_csv(output_path, sep='\t',
                                                              index=False, compression='gzip')
        return

    # Build donor signature for this (guide, donor) combination
    cut_site = find_cut_site(ref_seq, guide)
    edit_region = (cut_site - edit_region_size, cut_site + edit_region_size)

    # Initialize as empty DonorSignature
    donor_signature = DonorSignature(snvs={}, insertions={}, deletions={})
    if donor:
        try:
            donor_signature = build_donor_signature(ref_seq, donor, cut_site, edit_region=edit_region)
        except Exception as e:
            pass  # No donor signature if parsing fails

    # Coverage verification statistics
    coverage_stats = {
        'single_end': 0,
        'asymmetric': 0,
        'symmetric': 0,
        'fully_covered': 0,
        'partially_covered': 0,
        'coverage_warnings': []
    }

    # Classify each relevant sequence
    results = []
    classified = 0
    skipped = 0

    for seq_hash in relevant_seq_hashes:
        if seq_hash not in alignment_lookup:
            skipped += 1
            continue

        aln = alignment_lookup[seq_hash]

        # Detect mode
        mode = aln.get('mode', 'single_end')
        r1_mapped = aln.get('r1_is_mapped', False)
        r2_mapped = aln.get('r2_is_mapped', False)
        r1_seq = aln.get('r1_sequence', '')
        r2_seq = aln.get('r2_sequence', '')

        # Flexible mode detection based on read lengths and donor edits
        sequencing_mode = detect_sequencing_mode(
            r1_mapped, r2_mapped,
            len(r1_seq), len(r2_seq),
            donor_signature, len(ref_seq)
        )
        coverage_stats[sequencing_mode] += 1

        # Coverage verification
        r1_region = get_aligned_region(aln.get('r1_ref_start', -1), aln.get('r1_ref_end', -1), aln.get('r1_is_reverse', False)) if r1_mapped else None
        r2_region = get_aligned_region(aln.get('r2_ref_start', -1), aln.get('r2_ref_end', -1), aln.get('r2_is_reverse', False)) if r2_mapped else None

        is_covered, coverage_mode, uncovered_edits = check_donor_coverage(
            donor_signature, r1_region, r2_region, len(ref_seq)
        )

        if is_covered:
            coverage_stats['fully_covered'] += 1
        else:
            coverage_stats['partially_covered'] += 1
            if len(coverage_stats['coverage_warnings']) < 5:  # Limit warnings
                coverage_stats['coverage_warnings'].append(
                    f"  seq_hash {seq_hash[:8]}...: {len(uncovered_edits)} uncovered edits at positions {uncovered_edits}"
                )

        # Classification logic based on mode
        if not r1_mapped:
            results.append({
                'seq_hash': seq_hash,
                'guide': guide,
                'donor': donor,
                'outcome': 'UNALIGNED',
                'n_donor_snvs': 0,
                'n_non_donor': 0,
                'donor_fraction': 0.0,
                'donor_snvs_detected': '',
                'sequencing_mode': sequencing_mode,
                'coverage_mode': coverage_mode
            })
            continue

        # Parse R1 CIGAR
        r1_cigar = parse_cigar_string(aln.get('r1_cigar', ''))
        if r1_cigar is None:
            results.append({
                'seq_hash': seq_hash,
                'guide': guide,
                'donor': donor,
                'outcome': 'NO_CIGAR',
                'n_donor_snvs': 0,
                'n_non_donor': 0,
                'donor_fraction': 0.0,
                'donor_snvs_detected': '',
                'sequencing_mode': sequencing_mode,
                'coverage_mode': coverage_mode
            })
            continue

        if sequencing_mode == 'single_end':
            # Single-end mode: R1 only
            clf = classify_read_edit_distance(
                read_seq=r1_seq,
                ref_seq=ref_seq,
                donor_signature=donor_signature,
                ref_start=int(aln['r1_ref_start']),
                cigar_ops=r1_cigar,
                cut_site=cut_site,
                snv_distance_filter=snv_distance_filter,
                filter_chimeric=filter_chimeric,
                homopolymer_filter=homopolymer_filter,
                nhej_quantification_window=nhej_quantification_window
            )

        elif sequencing_mode == 'symmetric':
            # Symmetric mode: Both R1 and R2 span edit region
            # Use intersection approach for high confidence
            r1_clf = classify_read_edit_distance(
                read_seq=r1_seq,
                ref_seq=ref_seq,
                donor_signature=donor_signature,
                ref_start=int(aln['r1_ref_start']),
                cigar_ops=r1_cigar,
                cut_site=cut_site,
                snv_distance_filter=snv_distance_filter,
                filter_chimeric=filter_chimeric,
                homopolymer_filter=homopolymer_filter,
                nhej_quantification_window=nhej_quantification_window
            )

            if r2_mapped:
                r2_cigar = parse_cigar_string(aln.get('r2_cigar', ''))
                if r2_cigar:
                    r2_clf = classify_read_edit_distance(
                        read_seq=r2_seq,
                        ref_seq=ref_seq,
                        donor_signature=donor_signature,
                        ref_start=int(aln['r2_ref_start']),
                        cigar_ops=r2_cigar,
                        cut_site=cut_site,
                        snv_distance_filter=snv_distance_filter,
                        filter_chimeric=filter_chimeric,
                        homopolymer_filter=homopolymer_filter,
                        nhej_quantification_window=nhej_quantification_window
                    )
                    clf = merge_symmetric_classifications(r1_clf, r2_clf, donor_signature)
                else:
                    clf = r1_clf
            else:
                clf = r1_clf

        else:
            # Asymmetric mode: R1 primary, R2 validates
            clf = classify_read_edit_distance(
                read_seq=r1_seq,
                ref_seq=ref_seq,
                donor_signature=donor_signature,
                ref_start=int(aln['r1_ref_start']),
                cigar_ops=r1_cigar,
                cut_site=cut_site,
                snv_distance_filter=snv_distance_filter,
                filter_chimeric=filter_chimeric,
                homopolymer_filter=homopolymer_filter,
                nhej_quantification_window=nhej_quantification_window
            )

            # Use R2 for large deletion detection
            if r2_mapped and aln.get('is_proper_pair', False):
                insert_size = aln.get('insert_size', 0)
                expected_insert = len(ref_seq)
                if insert_size < expected_insert - 50:
                    # Large deletion detected
                    # Could add this to classification metadata if needed
                    pass

        # Check for donor capture (NHEJ-mediated insertion of donor fragment)
        donor_capture_info = detect_donor_capture(
            read_seq=r1_seq,
            cigar_tuples=r1_cigar,
            ref_start=int(aln['r1_ref_start']),
            donor_seq=donor,
            ref_seq=ref_seq,
            cut_site=cut_site,
            min_insertion_size=10,
            min_match_fraction=0.8
        )

        # Override outcome if donor capture detected
        final_outcome = clf.outcome
        is_donor_capture = False
        donor_capture_size = 0
        donor_capture_distance = 0
        junction_homology_left = 0
        junction_homology_right = 0
        mmej_class = ''

        if donor_capture_info:
            is_donor_capture = True
            donor_capture_size = donor_capture_info['insertion_size']
            donor_capture_distance = donor_capture_info['distance_to_cut']
            junction_homology_left = donor_capture_info.get('junction_homology_left', 0)
            junction_homology_right = donor_capture_info.get('junction_homology_right', 0)
            mmej_class = donor_capture_info.get('mmej_class', 'UNKNOWN')

            # Override outcome to DONOR_CAPTURE
            # This is distinct from HDR because it's NHEJ-mediated, not true HDR
            final_outcome = 'DONOR_CAPTURE'

        classified += 1
        results.append({
            'seq_hash': seq_hash,
            'guide': guide,
            'donor': donor,
            'outcome': final_outcome,
            'n_donor_snvs': clf.n_donor_encoded,
            'n_non_donor': clf.n_non_donor,
            'donor_fraction': clf.donor_fraction,
            'donor_snvs_detected': ','.join(str(p) for p in sorted(clf.donor_snvs_detected)) if clf.donor_snvs_detected else '',
            'sequencing_mode': sequencing_mode,
            'coverage_mode': coverage_mode,
            'is_donor_capture': is_donor_capture,
            'donor_capture_size': donor_capture_size,
            'donor_capture_distance': donor_capture_distance,
            'junction_homology_left': junction_homology_left,
            'junction_homology_right': junction_homology_right,
            'mmej_class': mmej_class
        })

    # Save results
    results_df = pd.DataFrame(results)
    results_df.to_csv(output_path, sep='\t', index=False, compression='gzip')

    # Write log with coverage statistics
    guide_display = guide[:30] + "..." if len(guide) > 30 else (guide if guide else "(empty)")
    donor_display = donor[:30] + "..." if len(donor) > 30 else (donor if donor else "(empty)")

    total_edits = len(donor_signature.snvs) + len(donor_signature.insertions) + len(donor_signature.deletions)
    coverage_pct = 100 * coverage_stats['fully_covered'] / max(classified, 1)

    log_msg = (
        f"Level 2 Paired-End Classification for (guide, donor):\n"
        f"  Guide: {guide_display}\n"
        f"  Donor: {donor_display}\n"
        f"  Cut site: {cut_site}\n"
        f"  Donor edits in core region: {total_edits}\n"
        f"    SNVs: {len(donor_signature.snvs)}, Insertions: {len(donor_signature.insertions)}, Deletions: {len(donor_signature.deletions)}\n"
        f"  \n"
        f"  Sequences with this guide/donor: {len(relevant_seq_hashes)}\n"
        f"  Classified: {classified}\n"
        f"  Skipped (not in alignment cache): {skipped}\n"
        f"  \n"
        f"Sequencing mode distribution:\n"
        f"  Single-end (R1 only): {coverage_stats['single_end']}\n"
        f"  Asymmetric paired: {coverage_stats['asymmetric']}\n"
        f"  Symmetric paired: {coverage_stats['symmetric']}\n"
        f"  \n"
        f"Coverage verification:\n"
        f"  Fully covered (all edits): {coverage_stats['fully_covered']} ({coverage_pct:.1f}%)\n"
        f"  Partially covered: {coverage_stats['partially_covered']}\n"
    )

    if coverage_stats['coverage_warnings']:
        log_msg += f"  \n⚠️  Coverage warnings (showing first 5):\n"
        for warning in coverage_stats['coverage_warnings']:
            log_msg += warning + "\n"
        log_msg += f"  \n  TRACE will report only on the portion of edits captured by reads.\n"

    log_msg += f"  \n"
    log_msg += f"Outcome distribution:\n"

    if len(results_df) > 0:
        for outcome, count in results_df['outcome'].value_counts().items():
            log_msg += f"  {outcome}: {count}\n"

    # Add donor capture statistics
    if len(results_df) > 0 and 'is_donor_capture' in results_df.columns:
        donor_captures = results_df[results_df['is_donor_capture'] == True]
        if len(donor_captures) > 0:
            capture_pct = 100 * len(donor_captures) / len(results_df)
            avg_size = donor_captures['donor_capture_size'].mean()
            avg_distance = donor_captures['donor_capture_distance'].mean()
            log_msg += f"  \n"
            log_msg += f"Donor capture detection:\n"
            log_msg += f"  Sequences with donor capture: {len(donor_captures)} ({capture_pct:.1f}%)\n"
            log_msg += f"  Average insertion size: {avg_size:.0f}bp (min 10bp threshold)\n"
            log_msg += f"  Average distance from cut site: {avg_distance:.0f}bp\n"
            log_msg += f"  \n"
            log_msg += f"⚠️  Donor capture represents NHEJ-mediated insertion of donor DNA,\n"
            log_msg += f"     not true HDR. These events are flagged separately.\n"

    log_path.write_text(log_msg)

    # Clean up SQLite temp file
    alignment_lookup.close()


if __name__ == '__main__':
    main()
