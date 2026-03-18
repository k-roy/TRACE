#!/usr/bin/env python3
"""
Classify all samples for a single (guide, donor) group.

Memory-efficient approach:
- Loads Level 1 alignment cache via SQLite (50MB vs 12GB)
- Builds Level 2 classification cache incrementally in memory
- Processes samples sequentially with checkpoints for resumability

Usage:
    python classify_group.py \
        --combo-id <hash> \
        --guide <sequence> \
        --donor <sequence> \
        --reference <sequence> \
        --manifest <path> \
        --alignment-cache <path> \
        --output-dir <path> \
        [--checkpoint-dir <path>]

Author: Kevin R. Roy
"""

import sys
import re
import argparse
import json
import tempfile
import sqlite3
import multiprocessing as mp
from functools import partial
from pathlib import Path
from collections import defaultdict
from dataclasses import dataclass
from typing import Optional, Set, Dict, List, Tuple

# Add TRACE to path
TRACE_ROOT = Path("/oak/stanford/groups/larsms/Users/kevinroy/software/trace")
sys.path.insert(0, str(TRACE_ROOT))

# Import fast alignment lookup (pre-indexed for O(1) lookups)
from alignment_cache_index import FastAlignmentLookup, AlignmentCacheIndex

import pandas as pd

from trace_crispr.config import parse_sequence_input
from trace_crispr.core.edit_distance_hdr import (
    build_donor_signature, classify_read_edit_distance, DonorSignature,
    classify_read_multi_donor, MultiDonorClassification, summarize_sample_swap_detection
)


def reverse_complement(seq):
    """Return reverse complement of DNA sequence."""
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
            'a': 't', 't': 'a', 'g': 'c', 'c': 'g', 'N': 'N', 'n': 'n'}
    return ''.join(comp.get(b, b) for b in reversed(seq))


def find_cut_site(ref_seq, guide_seq):
    """Find cut site position from guide location in reference."""
    if not guide_seq:
        return len(ref_seq) // 2

    ref_upper = ref_seq.upper()
    guide_upper = guide_seq.upper()

    # Try plus strand first
    guide_pos = ref_upper.find(guide_upper)
    if guide_pos != -1:
        return guide_pos + len(guide_seq) - 3

    # Try minus strand (reverse complement)
    guide_rc = reverse_complement(guide_upper)
    guide_pos = ref_upper.find(guide_rc)
    if guide_pos != -1:
        return guide_pos + 3

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
    """Analyze microhomology at insertion junctions to distinguish NHEJ vs MMEJ."""
    ins_upper = insertion_seq.upper()
    left_upper = ref_left_flank.upper()
    right_upper = ref_right_flank.upper()

    # Check left junction homology
    homology_left = 0
    for i in range(min(len(ins_upper), len(left_upper), max_homology_check)):
        if ins_upper[i] == left_upper[-(i+1)] if i < len(left_upper) else False:
            homology_left += 1
        else:
            break

    # Check right junction homology
    homology_right = 0
    for i in range(min(len(ins_upper), len(right_upper), max_homology_check)):
        if ins_upper[-(i+1)] == right_upper[i] if i < len(right_upper) else False:
            homology_right += 1
        else:
            break

    # MMEJ classification based on homology length
    max_homology = max(homology_left, homology_right)
    if max_homology >= 5:
        mmej_class = 'MMEJ_STRONG'
    elif max_homology >= 2:
        mmej_class = 'MMEJ_WEAK'
    else:
        mmej_class = 'NHEJ'

    return {
        'homology_left': homology_left,
        'homology_right': homology_right,
        'mmej_class': mmej_class
    }


def detect_donor_capture(read_seq, cigar_tuples, ref_start, donor_seq, ref_seq,
                         cut_site, min_insertion_size=10, min_match_fraction=0.8):
    """Detect if an insertion is a donor capture event."""
    if not cigar_tuples or not donor_seq:
        return None

    donor_upper = donor_seq.upper()
    insertions_near_cut = []
    ref_pos = ref_start
    read_pos = 0

    for op, length in cigar_tuples:
        if op == 0:  # M (match/mismatch)
            ref_pos += length
            read_pos += length
        elif op == 1:  # I (insertion)
            if length >= min_insertion_size:
                distance_to_cut = abs(ref_pos - cut_site)
                if distance_to_cut <= 50:
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

    for ins in insertions_near_cut:
        ins_seq = ins['seq'].upper()
        best_match_start = -1
        best_match_len = 0

        for start in range(len(donor_upper) - min_insertion_size + 1):
            for match_len in range(min(len(ins_seq), len(donor_upper) - start),
                                   min_insertion_size - 1, -1):
                donor_region = donor_upper[start:start + match_len]
                ins_region = ins_seq[:match_len]
                matches = sum(1 for a, b in zip(donor_region, ins_region) if a == b)
                match_frac = matches / match_len

                if match_frac >= min_match_fraction and match_len > best_match_len:
                    best_match_start = start
                    best_match_len = match_len

        if best_match_start >= 0 and best_match_len >= min_insertion_size:
            # Found donor capture - analyze junction microhomology
            ref_upper = ref_seq.upper() if ref_seq else ''
            ref_pos_in_seq = ins['ref_pos'] - ref_start
            flank_len = 15

            if ref_upper and 0 <= ref_pos_in_seq <= len(ref_upper):
                left_start = max(0, ref_pos_in_seq - flank_len)
                left_flank = ref_upper[left_start:ref_pos_in_seq]
                right_end = min(len(ref_upper), ref_pos_in_seq + flank_len)
                right_flank = ref_upper[ref_pos_in_seq:right_end]
            else:
                left_flank = ''
                right_flank = ''

            junction_info = analyze_junction_microhomology(
                insertion_seq=ins['seq'],
                ref_left_flank=left_flank,
                ref_right_flank=right_flank
            )

            return {
                'is_donor_capture': True,
                'insertion_size': ins['size'],
                'distance_to_cut': ins['distance_to_cut'],
                'junction_homology_left': junction_info['homology_left'],
                'junction_homology_right': junction_info['homology_right'],
                'mmej_class': junction_info['mmej_class']
            }

    return None


def get_aligned_region(ref_start, ref_end, is_reverse):
    """Get the reference region covered by a read alignment."""
    if ref_start < 0 or ref_end < 0:
        return None
    return (ref_start, ref_end)


def check_donor_coverage(donor_signature, r1_region, r2_region, ref_length):
    """Check if reads cover all donor edits."""
    edit_positions = set()
    edit_positions.update(donor_signature.snvs.keys())
    edit_positions.update(donor_signature.insertions.keys())
    edit_positions.update(donor_signature.deletions.keys())

    if not edit_positions:
        return True, 'no_edits', []

    r1_covered = set()
    if r1_region:
        r1_start, r1_end = r1_region
        r1_covered = {pos for pos in edit_positions if r1_start <= pos < r1_end}

    r2_covered = set()
    if r2_region:
        r2_start, r2_end = r2_region
        r2_covered = {pos for pos in edit_positions if r2_start <= pos < r2_end}

    all_covered = r1_covered | r2_covered
    uncovered = edit_positions - all_covered

    if len(uncovered) == 0:
        if len(r1_covered) == len(edit_positions):
            coverage_mode = 'r1_only'
        elif len(r2_covered) == len(edit_positions):
            coverage_mode = 'r2_only'
        else:
            coverage_mode = 'r1_and_r2'
        return True, coverage_mode, []
    else:
        return False, 'incomplete', sorted(uncovered)


def detect_sequencing_mode(r1_mapped, r2_mapped, r1_seq_len, r2_seq_len,
                           donor_signature, ref_length):
    """Detect sequencing mode based on read lengths."""
    if not r2_mapped:
        return 'single_end'

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

    if edit_span > 0:
        min_read_length = edit_span * 1.2
        if r1_seq_len >= min_read_length and r2_seq_len >= min_read_length:
            return 'symmetric'

    if r1_seq_len >= 200 and r2_seq_len >= 200:
        return 'symmetric'

    return 'asymmetric'


def merge_symmetric_classifications(r1_clf, r2_clf, donor_signature):
    """Merge R1 and R2 classifications for symmetric sequencing."""
    r1_donor_snvs = set(r1_clf.donor_snvs_detected) if r1_clf.donor_snvs_detected else set()
    r2_donor_snvs = set(r2_clf.donor_snvs_detected) if r2_clf.donor_snvs_detected else set()
    confirmed_donor_snvs = r1_donor_snvs & r2_donor_snvs

    if r1_clf.outcome == r2_clf.outcome:
        return r1_clf

    if r1_clf.outcome == 'WT':
        return r2_clf
    if r2_clf.outcome == 'WT':
        return r1_clf

    if 'HDR' in r1_clf.outcome and 'HDR' in r2_clf.outcome:
        total_donor_edits = len(donor_signature.snvs) + len(donor_signature.insertions) + len(donor_signature.deletions)
        if len(confirmed_donor_snvs) == total_donor_edits:
            r1_clf.outcome = 'HDR_COMPLETE'
            r1_clf.donor_snvs_detected = sorted(confirmed_donor_snvs)
            return r1_clf
        elif len(confirmed_donor_snvs) > 0:
            r1_clf.outcome = 'HDR_PARTIAL'
            r1_clf.donor_snvs_detected = sorted(confirmed_donor_snvs)
            return r1_clf

    return r1_clf


class SQLiteAlignmentLookup:
    """Memory-efficient alignment lookup using SQLite."""

    def __init__(self, alignment_cache_path: Path, relevant_hashes: Set[str]):
        """
        Create SQLite database from alignment cache TSV.

        Only loads sequences that are in relevant_hashes.
        Memory usage: ~50MB vs 4-12GB for full in-memory dict.
        """
        self.temp_db = tempfile.NamedTemporaryFile(suffix='.db', delete=False)
        self.temp_db.close()

        self.conn = sqlite3.connect(self.temp_db.name)
        self.conn.execute("PRAGMA journal_mode = WAL")
        self.conn.execute("PRAGMA synchronous = NORMAL")

        # Get columns
        first_chunk = pd.read_csv(
            alignment_cache_path,
            sep='\t',
            compression='gzip',
            nrows=1
        )
        all_columns = first_chunk.columns.tolist()
        self.columns = [c for c in all_columns if c != 'seq_hash']

        # Create table
        col_defs = ', '.join([f'"{col}" TEXT' for col in self.columns])
        self.conn.execute(f'CREATE TABLE alignments (seq_hash TEXT PRIMARY KEY, {col_defs})')

        # Load relevant hashes in chunks
        loaded = 0
        chunk_iter = pd.read_csv(
            alignment_cache_path,
            sep='\t',
            compression='gzip',
            chunksize=10000
        )

        for chunk in chunk_iter:
            relevant_chunk = chunk[chunk['seq_hash'].isin(relevant_hashes)]
            if len(relevant_chunk) > 0:
                insert_cols = ['seq_hash'] + self.columns
                placeholders = ', '.join(['?'] * len(insert_cols))

                batch_data = []
                for _, row in relevant_chunk.iterrows():
                    values = tuple(str(row[col]) if row[col] is not None else '' for col in insert_cols)
                    batch_data.append(values)

                self.conn.executemany(
                    f'INSERT OR REPLACE INTO alignments VALUES ({placeholders})',
                    batch_data
                )
                loaded += len(relevant_chunk)

        self.conn.commit()
        self.loaded_count = loaded

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

        col_names = [desc[0] for desc in cursor.description]
        result = {}
        for i, col in enumerate(col_names):
            val = row[i]
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
        Path(self.temp_db.name).unlink(missing_ok=True)


@dataclass
class ClassificationResult:
    """Container for a single sequence's classification."""
    seq_hash: str
    outcome: str
    n_donor_snvs: int
    n_non_donor: int
    donor_fraction: float
    donor_snvs_detected: str
    sequencing_mode: str
    coverage_mode: str
    uncovered_positions: str
    is_donor_capture: bool
    donor_capture_size: int
    donor_capture_distance: int
    junction_homology_left: int
    junction_homology_right: int
    mmej_class: str
    # Multi-donor comparison fields
    best_donor_id: str = ''
    expected_donor_id: str = ''
    is_sample_swap: bool = False
    swap_confidence: float = 0.0
    best_donor_score: float = 0.0


def classify_sequence_worker(args: Tuple) -> Tuple[str, ClassificationResult]:
    """
    Worker function for parallel classification.

    Takes a tuple of (seq_hash, aln_dict, shared_params) and returns (seq_hash, ClassificationResult).
    """
    seq_hash, aln, ref_seq, guide, donor, donor_sig_dict, cut_site, snv_distance_filter = args

    # Reconstruct DonorSignature from dict (needed for pickling across processes)
    donor_signature = DonorSignature(
        snvs=donor_sig_dict['snvs'],
        insertions=donor_sig_dict['insertions'],
        deletions=donor_sig_dict['deletions']
    )

    if aln is None:
        clf = ClassificationResult(
            seq_hash=seq_hash,
            outcome='NOT_IN_CACHE',
            n_donor_snvs=0,
            n_non_donor=0,
            donor_fraction=0.0,
            donor_snvs_detected='',
            sequencing_mode='unknown',
            coverage_mode='unknown',
            uncovered_positions='',
            is_donor_capture=False,
            donor_capture_size=0,
            donor_capture_distance=0,
            junction_homology_left=0,
            junction_homology_right=0,
            mmej_class=''
        )
    else:
        clf = classify_sequence(
            seq_hash=seq_hash,
            aln=aln,
            ref_seq=ref_seq,
            guide=guide,
            donor=donor,
            donor_signature=donor_signature,
            cut_site=cut_site,
            snv_distance_filter=snv_distance_filter
        )

    return (seq_hash, clf)


def classify_sequence_multi_donor_worker(args: Tuple) -> Tuple[str, ClassificationResult]:
    """
    Worker function for parallel multi-donor classification.

    Compares each sequence against all donor templates to find the best match
    and detect potential sample swaps.
    """
    (seq_hash, aln, ref_seq, guide, expected_donor_id,
     all_donor_sig_dicts, cut_site, snv_distance_filter) = args

    # Reconstruct all DonorSignatures from dicts
    all_donor_signatures = {}
    for donor_id, sig_dict in all_donor_sig_dicts.items():
        all_donor_signatures[donor_id] = DonorSignature(
            snvs=sig_dict['snvs'],
            insertions=sig_dict['insertions'],
            deletions=sig_dict['deletions']
        )

    if aln is None:
        return (seq_hash, ClassificationResult(
            seq_hash=seq_hash,
            outcome='NOT_IN_CACHE',
            n_donor_snvs=0,
            n_non_donor=0,
            donor_fraction=0.0,
            donor_snvs_detected='',
            sequencing_mode='unknown',
            coverage_mode='unknown',
            uncovered_positions='',
            is_donor_capture=False,
            donor_capture_size=0,
            donor_capture_distance=0,
            junction_homology_left=0,
            junction_homology_right=0,
            mmej_class='',
            best_donor_id='',
            expected_donor_id=expected_donor_id,
            is_sample_swap=False,
            swap_confidence=0.0,
            best_donor_score=0.0
        ))

    clf = classify_sequence_multi_donor(
        seq_hash=seq_hash,
        aln=aln,
        ref_seq=ref_seq,
        guide=guide,
        expected_donor_id=expected_donor_id,
        all_donor_signatures=all_donor_signatures,
        cut_site=cut_site,
        snv_distance_filter=snv_distance_filter
    )

    return (seq_hash, clf)


def classify_sequence_multi_donor(
    seq_hash: str, aln: dict, ref_seq: str, guide: str,
    expected_donor_id: str, all_donor_signatures: Dict[str, DonorSignature],
    cut_site: int, snv_distance_filter: int = 50
) -> ClassificationResult:
    """
    Classify a single sequence against multiple donor templates.

    Returns ClassificationResult with best matching donor and swap detection.
    """
    mode = aln.get('mode', 'single_end')
    r1_mapped = aln.get('r1_is_mapped', False)
    r2_mapped = aln.get('r2_is_mapped', False)
    r1_seq = aln.get('r1_sequence', '')
    r2_seq = aln.get('r2_sequence', '')

    # Use expected donor for sequencing mode detection (representative)
    expected_signature = all_donor_signatures.get(expected_donor_id, DonorSignature({}, {}, {}))

    sequencing_mode = detect_sequencing_mode(
        r1_mapped, r2_mapped,
        len(r1_seq), len(r2_seq),
        expected_signature, len(ref_seq)
    )

    # Coverage verification using expected donor
    r1_region = get_aligned_region(
        aln.get('r1_ref_start', -1), aln.get('r1_ref_end', -1),
        aln.get('r1_is_reverse', False)
    ) if r1_mapped else None
    r2_region = get_aligned_region(
        aln.get('r2_ref_start', -1), aln.get('r2_ref_end', -1),
        aln.get('r2_is_reverse', False)
    ) if r2_mapped else None

    is_covered, coverage_mode, uncovered_edits = check_donor_coverage(
        expected_signature, r1_region, r2_region, len(ref_seq)
    )
    uncovered_str = ','.join(str(p) for p in uncovered_edits) if uncovered_edits else ''

    # Handle unmapped reads
    if not r1_mapped:
        return ClassificationResult(
            seq_hash=seq_hash,
            outcome='UNALIGNED',
            n_donor_snvs=0,
            n_non_donor=0,
            donor_fraction=0.0,
            donor_snvs_detected='',
            sequencing_mode=sequencing_mode,
            coverage_mode=coverage_mode,
            uncovered_positions=uncovered_str,
            is_donor_capture=False,
            donor_capture_size=0,
            donor_capture_distance=0,
            junction_homology_left=0,
            junction_homology_right=0,
            mmej_class='',
            best_donor_id='',
            expected_donor_id=expected_donor_id,
            is_sample_swap=False,
            swap_confidence=0.0,
            best_donor_score=0.0
        )

    # Parse CIGAR
    r1_cigar = parse_cigar_string(aln.get('r1_cigar', ''))
    if r1_cigar is None:
        return ClassificationResult(
            seq_hash=seq_hash,
            outcome='NO_CIGAR',
            n_donor_snvs=0,
            n_non_donor=0,
            donor_fraction=0.0,
            donor_snvs_detected='',
            sequencing_mode=sequencing_mode,
            coverage_mode=coverage_mode,
            uncovered_positions=uncovered_str,
            is_donor_capture=False,
            donor_capture_size=0,
            donor_capture_distance=0,
            junction_homology_left=0,
            junction_homology_right=0,
            mmej_class='',
            best_donor_id='',
            expected_donor_id=expected_donor_id,
            is_sample_swap=False,
            swap_confidence=0.0,
            best_donor_score=0.0
        )

    # Multi-donor classification
    multi_clf = classify_read_multi_donor(
        read_seq=r1_seq,
        ref_seq=ref_seq,
        donor_signatures=all_donor_signatures,
        expected_donor_id=expected_donor_id,
        ref_start=int(aln['r1_ref_start']),
        cigar_ops=r1_cigar,
        cut_site=cut_site,
        snv_distance_filter=snv_distance_filter
    )

    # Get classification from best matching donor
    clf = multi_clf.best_donor_classification

    # Check for donor capture using expected donor
    donor_capture_info = detect_donor_capture(
        read_seq=r1_seq,
        cigar_tuples=r1_cigar,
        ref_start=int(aln['r1_ref_start']),
        donor_seq='',  # Would need to reconstruct from signature
        ref_seq=ref_seq,
        cut_site=cut_site,
        min_insertion_size=10,
        min_match_fraction=0.8
    )

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
        final_outcome = 'DONOR_CAPTURE'

    donor_snvs_str = ','.join(str(p) for p in sorted(clf.donor_snvs_detected)) if clf.donor_snvs_detected else ''

    return ClassificationResult(
        seq_hash=seq_hash,
        outcome=final_outcome,
        n_donor_snvs=clf.n_donor_encoded,
        n_non_donor=clf.n_non_donor,
        donor_fraction=clf.donor_fraction,
        donor_snvs_detected=donor_snvs_str,
        sequencing_mode=sequencing_mode,
        coverage_mode=coverage_mode,
        uncovered_positions=uncovered_str,
        is_donor_capture=is_donor_capture,
        donor_capture_size=donor_capture_size,
        donor_capture_distance=donor_capture_distance,
        junction_homology_left=junction_homology_left,
        junction_homology_right=junction_homology_right,
        mmej_class=mmej_class,
        best_donor_id=multi_clf.best_donor_id,
        expected_donor_id=multi_clf.expected_donor_id,
        is_sample_swap=multi_clf.is_sample_swap,
        swap_confidence=multi_clf.swap_confidence,
        best_donor_score=multi_clf.best_donor_score
    )


def classify_sequence(seq_hash: str, aln: dict, ref_seq: str, guide: str, donor: str,
                      donor_signature: DonorSignature, cut_site: int,
                      snv_distance_filter: int = 50) -> ClassificationResult:
    """
    Classify a single sequence using cached alignment.

    Returns ClassificationResult with outcome and metadata.
    """
    mode = aln.get('mode', 'single_end')
    r1_mapped = aln.get('r1_is_mapped', False)
    r2_mapped = aln.get('r2_is_mapped', False)
    r1_seq = aln.get('r1_sequence', '')
    r2_seq = aln.get('r2_sequence', '')

    # Detect sequencing mode
    sequencing_mode = detect_sequencing_mode(
        r1_mapped, r2_mapped,
        len(r1_seq), len(r2_seq),
        donor_signature, len(ref_seq)
    )

    # Coverage verification
    r1_region = get_aligned_region(
        aln.get('r1_ref_start', -1), aln.get('r1_ref_end', -1),
        aln.get('r1_is_reverse', False)
    ) if r1_mapped else None
    r2_region = get_aligned_region(
        aln.get('r2_ref_start', -1), aln.get('r2_ref_end', -1),
        aln.get('r2_is_reverse', False)
    ) if r2_mapped else None

    is_covered, coverage_mode, uncovered_edits = check_donor_coverage(
        donor_signature, r1_region, r2_region, len(ref_seq)
    )
    uncovered_str = ','.join(str(p) for p in uncovered_edits) if uncovered_edits else ''

    # Handle unmapped reads
    if not r1_mapped:
        return ClassificationResult(
            seq_hash=seq_hash,
            outcome='UNALIGNED',
            n_donor_snvs=0,
            n_non_donor=0,
            donor_fraction=0.0,
            donor_snvs_detected='',
            sequencing_mode=sequencing_mode,
            coverage_mode=coverage_mode,
            uncovered_positions=uncovered_str,
            is_donor_capture=False,
            donor_capture_size=0,
            donor_capture_distance=0,
            junction_homology_left=0,
            junction_homology_right=0,
            mmej_class=''
        )

    # Parse CIGAR
    r1_cigar = parse_cigar_string(aln.get('r1_cigar', ''))
    if r1_cigar is None:
        return ClassificationResult(
            seq_hash=seq_hash,
            outcome='NO_CIGAR',
            n_donor_snvs=0,
            n_non_donor=0,
            donor_fraction=0.0,
            donor_snvs_detected='',
            sequencing_mode=sequencing_mode,
            coverage_mode=coverage_mode,
            uncovered_positions=uncovered_str,
            is_donor_capture=False,
            donor_capture_size=0,
            donor_capture_distance=0,
            junction_homology_left=0,
            junction_homology_right=0,
            mmej_class=''
        )

    # Classification based on mode
    if sequencing_mode == 'single_end':
        clf = classify_read_edit_distance(
            read_seq=r1_seq,
            ref_seq=ref_seq,
            donor_signature=donor_signature,
            ref_start=int(aln['r1_ref_start']),
            cigar_ops=r1_cigar,
            cut_site=cut_site,
            snv_distance_filter=snv_distance_filter
        )
    elif sequencing_mode == 'symmetric':
        r1_clf = classify_read_edit_distance(
            read_seq=r1_seq,
            ref_seq=ref_seq,
            donor_signature=donor_signature,
            ref_start=int(aln['r1_ref_start']),
            cigar_ops=r1_cigar,
            cut_site=cut_site,
            snv_distance_filter=snv_distance_filter
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
                    snv_distance_filter=snv_distance_filter
                )
                clf = merge_symmetric_classifications(r1_clf, r2_clf, donor_signature)
            else:
                clf = r1_clf
        else:
            clf = r1_clf
    else:
        # Asymmetric mode
        clf = classify_read_edit_distance(
            read_seq=r1_seq,
            ref_seq=ref_seq,
            donor_signature=donor_signature,
            ref_start=int(aln['r1_ref_start']),
            cigar_ops=r1_cigar,
            cut_site=cut_site,
            snv_distance_filter=snv_distance_filter
        )

    # Check for donor capture
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
        final_outcome = 'DONOR_CAPTURE'

    donor_snvs_str = ','.join(str(p) for p in sorted(clf.donor_snvs_detected)) if clf.donor_snvs_detected else ''

    return ClassificationResult(
        seq_hash=seq_hash,
        outcome=final_outcome,
        n_donor_snvs=clf.n_donor_encoded,
        n_non_donor=clf.n_non_donor,
        donor_fraction=clf.donor_fraction,
        donor_snvs_detected=donor_snvs_str,
        sequencing_mode=sequencing_mode,
        coverage_mode=coverage_mode,
        uncovered_positions=uncovered_str,
        is_donor_capture=is_donor_capture,
        donor_capture_size=donor_capture_size,
        donor_capture_distance=donor_capture_distance,
        junction_homology_left=junction_homology_left,
        junction_homology_right=junction_homology_right,
        mmej_class=mmej_class
    )


def write_sample_outputs(sample_id: str, results: list, output_dir: Path,
                         multi_donor_mode: bool = False):
    """
    Write per-sample classification.tsv and hdr_snv_detail.tsv.

    Uses atomic writes (write to temp, then rename) for crash safety.
    """
    sample_dir = output_dir / "samples" / sample_id
    sample_dir.mkdir(parents=True, exist_ok=True)

    # Prepare results DataFrame
    rows = []
    per_snv_counts = defaultdict(int)
    total_hdr_reads = 0

    for seq_hash, clf, count in results:
        row = {
            'read_name': f"{seq_hash}_x{count}",
            'outcome': clf.outcome,
            'n_donor_snvs': clf.n_donor_snvs,
            'n_non_donor': clf.n_non_donor,
            'donor_fraction': clf.donor_fraction,
            'donor_snvs_detected': clf.donor_snvs_detected,
            'sequencing_mode': clf.sequencing_mode,
            'coverage_mode': clf.coverage_mode,
            'uncovered_positions': clf.uncovered_positions,
            'is_donor_capture': clf.is_donor_capture,
            'donor_capture_size': clf.donor_capture_size,
            'donor_capture_distance': clf.donor_capture_distance,
            'junction_homology_left': clf.junction_homology_left,
            'junction_homology_right': clf.junction_homology_right,
            'mmej_class': clf.mmej_class,
            'count': count
        }

        # Add multi-donor columns if enabled
        if multi_donor_mode:
            row['best_donor_id'] = clf.best_donor_id
            row['expected_donor_id'] = clf.expected_donor_id
            row['is_sample_swap'] = clf.is_sample_swap
            row['swap_confidence'] = clf.swap_confidence
            row['best_donor_score'] = clf.best_donor_score

        rows.append(row)

        # Track per-SNV integration for HDR reads
        if clf.outcome in ['HDR_COMPLETE', 'HDR_PARTIAL', 'MIXED']:
            total_hdr_reads += count
            if clf.donor_snvs_detected:
                for pos in clf.donor_snvs_detected.split(','):
                    if pos:
                        per_snv_counts[int(pos)] += count

    results_df = pd.DataFrame(rows)

    # Write classification.tsv atomically
    classification_path = sample_dir / "classification.tsv"
    temp_clf = classification_path.with_suffix('.tsv.tmp')
    results_df.to_csv(temp_clf, sep='\t', index=False)
    temp_clf.rename(classification_path)

    # Write HDR detail
    hdr_detail_rows = []
    for pos, count in sorted(per_snv_counts.items()):
        freq = count / total_hdr_reads if total_hdr_reads > 0 else 0
        hdr_detail_rows.append({
            'position': pos,
            'count': count,
            'frequency': freq
        })

    hdr_detail_df = pd.DataFrame(hdr_detail_rows)
    hdr_detail_path = sample_dir / "hdr_snv_detail.tsv"
    temp_hdr = hdr_detail_path.with_suffix('.tsv.tmp')
    hdr_detail_df.to_csv(temp_hdr, sep='\t', index=False)
    temp_hdr.rename(hdr_detail_path)

    return len(results), results_df['count'].sum() if len(rows) > 0 else 0


def load_checkpoint(checkpoint_path: Path) -> set:
    """Load set of completed sample IDs from checkpoint file."""
    if not checkpoint_path.exists():
        return set()

    completed = set()
    with open(checkpoint_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line:
                completed.add(line)
    return completed


def append_checkpoint(checkpoint_path: Path, sample_id: str):
    """Append completed sample ID to checkpoint file."""
    with open(checkpoint_path, 'a') as f:
        f.write(f"{sample_id}\n")


def main():
    parser = argparse.ArgumentParser(description='Classify all samples for a (guide, donor) group')
    parser.add_argument('--combo-id', required=True, help='Combo ID hash')
    parser.add_argument('--guide', default='', help='Guide sequence')
    parser.add_argument('--donor', default='', help='Donor/HDR template sequence')
    parser.add_argument('--reference', required=True, help='Reference sequence')
    parser.add_argument('--manifest', required=True, help='Path to sample manifest')
    parser.add_argument('--alignment-cache', required=True, help='Path to Level 1 alignment cache')
    parser.add_argument('--output-dir', required=True, help='Output directory')
    parser.add_argument('--samples-dir', required=True, help='Directory containing sample collapsed files')
    parser.add_argument('--checkpoint-dir', help='Checkpoint directory (default: output-dir/checkpoints)')
    parser.add_argument('--snv-distance-filter', type=int, default=50)
    parser.add_argument('--multi-donor-mode', action='store_true',
                        help='Enable multi-donor comparison for barcoded experiments. '
                             'Compares each read against all donor templates in the manifest '
                             'to find the best match and detect sample swaps.')

    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    samples_dir = Path(args.samples_dir)
    checkpoint_dir = Path(args.checkpoint_dir) if args.checkpoint_dir else output_dir / "checkpoints"
    checkpoint_dir.mkdir(parents=True, exist_ok=True)

    checkpoint_path = checkpoint_dir / f"{args.combo_id}.checkpoint"
    log_path = output_dir / "logs" / f"classify_group_{args.combo_id}.log"
    log_path.parent.mkdir(parents=True, exist_ok=True)

    # Parse reference
    ref_seq = parse_sequence_input(args.reference)

    # Parse guide and donor (may be empty)
    guide = args.guide if args.guide and args.guide != 'nan' else ''
    donor = args.donor if args.donor and args.donor != 'nan' else ''

    # Find cut site first (needed for donor signature)
    cut_site = find_cut_site(ref_seq, guide)

    # Build donor signature
    donor_signature = build_donor_signature(ref_seq, donor, cut_site) if donor else DonorSignature({}, {}, {})

    # Load manifest and find samples for this combo
    manifest_df = pd.read_csv(args.manifest, sep='\t')
    samples_for_combo = []
    sample_to_expected_donor = {}  # sample_id -> expected donor ID (for multi-donor mode)

    for _, row in manifest_df.iterrows():
        sample_guide = str(row.get('guide', '')) if pd.notna(row.get('guide', '')) else ''
        sample_donor = str(row.get('hdr_template', '')) if pd.notna(row.get('hdr_template', '')) else ''

        # Match guide/donor
        if sample_guide == guide and sample_donor == donor:
            samples_for_combo.append(row['sample_id'])
            sample_to_expected_donor[row['sample_id']] = sample_donor

    # Multi-donor mode: load all unique donors for this guide
    all_donor_signatures = {}  # donor_seq -> DonorSignature
    all_donor_ids = {}  # donor_seq -> short ID for reporting
    if args.multi_donor_mode:
        # Find all unique donors that share this guide
        unique_donors = set()
        for _, row in manifest_df.iterrows():
            sample_guide = str(row.get('guide', '')) if pd.notna(row.get('guide', '')) else ''
            sample_donor = str(row.get('hdr_template', '')) if pd.notna(row.get('hdr_template', '')) else ''

            if sample_guide == guide and sample_donor:
                unique_donors.add(sample_donor)

        print(f"Multi-donor mode: building signatures for {len(unique_donors)} unique donors...")
        for i, donor_seq in enumerate(sorted(unique_donors)):
            # Create a short ID for this donor (using first 8 chars of hash)
            import hashlib
            donor_hash = hashlib.md5(donor_seq.encode()).hexdigest()[:8]
            donor_id = f"donor_{donor_hash}"
            all_donor_ids[donor_seq] = donor_id

            # Build signature for this donor
            sig = build_donor_signature(ref_seq, donor_seq, cut_site)
            all_donor_signatures[donor_id] = sig

            if (i + 1) % 10 == 0 or i + 1 == len(unique_donors):
                print(f"  Built {i+1}/{len(unique_donors)} donor signatures")

        # Map sample expected donors to IDs
        for sample_id in samples_for_combo:
            expected_donor_seq = sample_to_expected_donor.get(sample_id, '')
            if expected_donor_seq:
                sample_to_expected_donor[sample_id] = all_donor_ids.get(expected_donor_seq, '')
            else:
                sample_to_expected_donor[sample_id] = ''

    print(f"Found {len(samples_for_combo)} samples for combo {args.combo_id}")
    print(f"  Guide: {guide[:30]}..." if len(guide) > 30 else f"  Guide: {guide or '(none)'}")
    print(f"  Donor: {donor[:30]}..." if len(donor) > 30 else f"  Donor: {donor or '(none)'}")

    # Load checkpoint
    completed_samples = load_checkpoint(checkpoint_path)
    remaining_samples = [s for s in samples_for_combo if s not in completed_samples]
    print(f"  Completed: {len(completed_samples)}, Remaining: {len(remaining_samples)}")

    if not remaining_samples:
        print("All samples already completed!")
        return

    # Collect all unique seq_hashes across remaining samples
    print("Collecting unique sequences across samples...")
    all_seq_hashes = set()
    sample_sequences = {}  # sample_id -> list of (seq_hash, count)

    for sample_id in remaining_samples:
        collapsed_path = samples_dir / sample_id / "collapsed" / "unique_sequences.tsv.gz"
        if not collapsed_path.exists():
            print(f"  Warning: {collapsed_path} not found, skipping")
            continue

        sample_df = pd.read_csv(collapsed_path, sep='\t', compression='gzip')
        seqs = [(row['seq_hash'], row['count']) for _, row in sample_df.iterrows()]
        sample_sequences[sample_id] = seqs
        all_seq_hashes.update(h for h, _ in seqs)

    print(f"  Total unique sequences: {len(all_seq_hashes)}")

    # Load Level 1 alignment cache via pre-indexed lookup (O(1) per sequence)
    # First run builds index (~6 min), subsequent runs use cached index (~0s)
    print("Loading alignment cache via pre-indexed lookup...")
    cache_dir = Path(args.alignment_cache).parent
    alignment_lookup = FastAlignmentLookup(
        Path(args.alignment_cache),
        all_seq_hashes,
        cache_dir=cache_dir
    )

    # Determine number of parallel workers
    n_workers = min(mp.cpu_count(), 16)  # Cap at 16 to match typical SLURM allocation
    print(f"Using {n_workers} parallel workers for classification")

    # Prepare donor signature(s) as dict for pickling
    donor_sig_dict = {
        'snvs': dict(donor_signature.snvs),
        'insertions': dict(donor_signature.insertions),
        'deletions': dict(donor_signature.deletions)
    }

    # For multi-donor mode, prepare all donor signatures as dicts
    all_donor_sig_dicts = {}
    if args.multi_donor_mode:
        for donor_id, sig in all_donor_signatures.items():
            all_donor_sig_dicts[donor_id] = {
                'snvs': dict(sig.snvs),
                'insertions': dict(sig.insertions),
                'deletions': dict(sig.deletions)
            }

    # Extract alignments in batches (64K lookups/sec with pre-indexed cache)
    print("Extracting alignments for parallel classification...")
    seq_list = list(all_seq_hashes)
    seq_to_aln = {}

    # Batch fetch for efficiency
    batch_size = 10000
    for i in range(0, len(seq_list), batch_size):
        batch = seq_list[i:i+batch_size]
        batch_results = alignment_lookup.get_batch(batch)
        seq_to_aln.update(batch_results)
        if (i + batch_size) % 100000 == 0 or i + batch_size >= len(seq_list):
            print(f"  Fetched {min(i + batch_size, len(seq_list)):,}/{len(seq_list):,} alignments...")

    # Close lookup
    alignment_lookup.close()

    # Prepare work items for parallel classification
    if args.multi_donor_mode:
        # Multi-donor mode: classify against all donors (use '' as placeholder expected_donor)
        work_items = [
            (seq_hash, seq_to_aln.get(seq_hash), ref_seq, guide, '',
             all_donor_sig_dicts, cut_site, args.snv_distance_filter)
            for seq_hash in all_seq_hashes
        ]
        worker_func = classify_sequence_multi_donor_worker
    else:
        # Standard mode: single donor
        work_items = [
            (seq_hash, seq_to_aln.get(seq_hash), ref_seq, guide, donor, donor_sig_dict,
             cut_site, args.snv_distance_filter)
            for seq_hash in all_seq_hashes
        ]
        worker_func = classify_sequence_worker

    # Classify all unique sequences in parallel
    mode_str = "multi-donor" if args.multi_donor_mode else "single-donor"
    print(f"Classifying {len(work_items)} unique sequences in parallel ({mode_str} mode)...")
    level2_cache = {}  # type: Dict[str, ClassificationResult]

    with mp.Pool(processes=n_workers) as pool:
        # Use imap_unordered for better progress tracking
        results_iter = pool.imap_unordered(worker_func, work_items, chunksize=100)

        classified = 0
        for seq_hash, clf in results_iter:
            level2_cache[seq_hash] = clf
            classified += 1
            if classified % 10000 == 0:
                print(f"  Classified {classified}/{len(work_items)} sequences...")

    print(f"  Classification complete: {len(level2_cache)} sequences")

    # Free memory from alignment dict (no longer needed)
    del seq_to_aln

    # Now expand to samples (fast - just cache lookups and writes)
    print("Expanding classifications to samples...")
    log_lines = []
    total_classified = 0
    total_reads = 0
    swap_summary = []  # Track sample swaps for reporting

    for i, sample_id in enumerate(remaining_samples):
        if sample_id not in sample_sequences:
            continue

        seqs = sample_sequences[sample_id]
        results = []

        # Get expected donor for this sample (multi-donor mode)
        expected_donor_id = sample_to_expected_donor.get(sample_id, '') if args.multi_donor_mode else ''

        for seq_hash, count in seqs:
            clf = level2_cache.get(seq_hash)
            if clf is None:
                # Should not happen, but handle gracefully
                clf = ClassificationResult(
                    seq_hash=seq_hash,
                    outcome='MISSING',
                    n_donor_snvs=0,
                    n_non_donor=0,
                    donor_fraction=0.0,
                    donor_snvs_detected='',
                    sequencing_mode='unknown',
                    coverage_mode='unknown',
                    uncovered_positions='',
                    is_donor_capture=False,
                    donor_capture_size=0,
                    donor_capture_distance=0,
                    junction_homology_left=0,
                    junction_homology_right=0,
                    mmej_class=''
                )
            elif args.multi_donor_mode and expected_donor_id:
                # Update classification with sample-specific expected donor and swap detection
                is_swap = (clf.best_donor_id != expected_donor_id and
                           clf.best_donor_score > 0 and
                           clf.n_donor_snvs > 0)
                swap_confidence = clf.swap_confidence if is_swap else 0.0

                # Create updated classification with correct expected_donor_id
                clf = ClassificationResult(
                    seq_hash=clf.seq_hash,
                    outcome=clf.outcome,
                    n_donor_snvs=clf.n_donor_snvs,
                    n_non_donor=clf.n_non_donor,
                    donor_fraction=clf.donor_fraction,
                    donor_snvs_detected=clf.donor_snvs_detected,
                    sequencing_mode=clf.sequencing_mode,
                    coverage_mode=clf.coverage_mode,
                    uncovered_positions=clf.uncovered_positions,
                    is_donor_capture=clf.is_donor_capture,
                    donor_capture_size=clf.donor_capture_size,
                    donor_capture_distance=clf.donor_capture_distance,
                    junction_homology_left=clf.junction_homology_left,
                    junction_homology_right=clf.junction_homology_right,
                    mmej_class=clf.mmej_class,
                    best_donor_id=clf.best_donor_id,
                    expected_donor_id=expected_donor_id,
                    is_sample_swap=is_swap,
                    swap_confidence=swap_confidence,
                    best_donor_score=clf.best_donor_score
                )
            results.append((seq_hash, clf, count))

        # Write sample outputs
        n_seqs, n_reads = write_sample_outputs(sample_id, results, output_dir, args.multi_donor_mode)

        # Summarize sample swap detection for multi-donor mode
        if args.multi_donor_mode:
            swap_reads = sum(count for _, clf, count in results if clf.is_sample_swap)
            hdr_reads = sum(count for _, clf, count in results if clf.n_donor_snvs > 0)
            if hdr_reads > 0 and swap_reads > hdr_reads * 0.5:
                # Majority of HDR reads show different donor - likely swap
                best_donors = {}
                for _, clf, count in results:
                    if clf.n_donor_snvs > 0:
                        best_donors[clf.best_donor_id] = best_donors.get(clf.best_donor_id, 0) + count
                if best_donors:
                    detected_donor = max(best_donors, key=best_donors.get)
                    swap_summary.append({
                        'sample_id': sample_id,
                        'expected_donor': expected_donor_id,
                        'detected_donor': detected_donor,
                        'swap_reads': swap_reads,
                        'total_hdr_reads': hdr_reads,
                        'swap_fraction': swap_reads / hdr_reads
                    })
        total_classified += n_seqs
        total_reads += n_reads

        # Update checkpoint
        append_checkpoint(checkpoint_path, sample_id)

        # Progress logging
        progress = f"[{i+1}/{len(remaining_samples)}]"
        log_line = f"{progress} {sample_id}: {n_seqs} seqs, {n_reads} reads"
        print(log_line)
        log_lines.append(log_line)

    # Write swap detection report if in multi-donor mode
    swap_report_lines = []
    if args.multi_donor_mode and swap_summary:
        print(f"\n=== SAMPLE SWAP DETECTION ===")
        print(f"Found {len(swap_summary)} potential sample swaps:")
        swap_report_lines.append("sample_id\texpected_donor\tdetected_donor\tswap_reads\ttotal_hdr_reads\tswap_fraction")
        for swap in swap_summary:
            line = f"{swap['sample_id']}\t{swap['expected_donor']}\t{swap['detected_donor']}\t{swap['swap_reads']}\t{swap['total_hdr_reads']}\t{swap['swap_fraction']:.3f}"
            swap_report_lines.append(line)
            print(f"  {swap['sample_id']}: expected={swap['expected_donor']}, detected={swap['detected_donor']} ({swap['swap_fraction']:.1%})")

        # Write swap report
        swap_report_path = output_dir / "swap_detection_report.tsv"
        swap_report_path.write_text('\n'.join(swap_report_lines))
        print(f"Swap report written to: {swap_report_path}")

    # Write final log
    mode_str = "Multi-donor mode" if args.multi_donor_mode else "Single-donor mode"
    donor_info = f"Donors: {len(all_donor_signatures)}" if args.multi_donor_mode else f"Donor: {donor[:50] if donor else '(none)'}"

    log_msg = f"""Classify Group Complete: {args.combo_id}

Mode: {mode_str}
Guide: {guide[:50] if guide else '(none)'}
{donor_info}

Samples processed: {len(remaining_samples)}
Level 2 cache entries: {len(level2_cache)}
Total sequences classified: {total_classified}
Total reads: {total_reads}
"""

    if args.multi_donor_mode:
        log_msg += f"\nSample swaps detected: {len(swap_summary)}\n"
        if swap_summary:
            log_msg += "\nSwap details:\n"
            for swap in swap_summary:
                log_msg += f"  {swap['sample_id']}: {swap['expected_donor']} -> {swap['detected_donor']} ({swap['swap_fraction']:.1%})\n"

    log_msg += "\nSample details:\n" + '\n'.join(log_lines)

    log_path.write_text(log_msg)
    print(f"\nDone! Log written to {log_path}")


if __name__ == '__main__':
    main()
