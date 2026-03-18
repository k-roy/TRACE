#!/usr/bin/env python3
"""
Level 1 Cache: Extract globally unique R1/R2 PAIRS (ignoring guide/donor).

SQLite-based implementation for low memory usage (~50MB vs 22GB).

NEW APPROACH: Instead of merged sequences, track unique R1/R2 pairs.
This creates a deduplicated set of paired reads for alignment.

The key insight: alignment only depends on (r1_seq, r2_seq, reference),
NOT on guide or donor. So we align each unique R1/R2 pair once.

Supports three modes:
- Paired-end: Both R1 and R2 sequences present
- Single-end: Only R1 sequences (R2 is empty string)

Output:
- level1_unique_sequences.tsv.gz: unique R1/R2 pairs to align
- level1_seq_to_samples.tsv.gz: maps each pair to samples that contain it
- level1_cache.db: SQLite database for downstream processing

Memory usage: ~50MB regardless of dataset size (was 22GB+ for large datasets).
"""

import gzip
import sys
from pathlib import Path

import pandas as pd

# Add parent directory for imports
WORKFLOW_DIR = Path(__file__).parent.parent
sys.path.insert(0, str(WORKFLOW_DIR.parent))

from trace_crispr.io.sqlite_cache import SequenceCache


def process_collapsed_file(cache: SequenceCache, collapsed_file: Path, sample_id: str):
    """
    Process a collapsed sequences file and add to SQLite cache.

    Uses chunked reading to minimize memory usage.
    """
    # Read in chunks to avoid loading entire file into memory
    chunk_iter = pd.read_csv(
        collapsed_file,
        sep='\t',
        compression='gzip',
        chunksize=10000
    )

    for chunk in chunk_iter:
        # Convert chunk to iterator of dicts
        def row_iterator():
            for _, row in chunk.iterrows():
                yield {
                    'r1_sequence': row['r1_sequence'],
                    'r2_sequence': row['r2_sequence'] if pd.notna(row['r2_sequence']) else '',
                    'count': row['count'],
                    'guide': row['guide'] if pd.notna(row['guide']) else '',
                    'donor': row['donor'] if pd.notna(row['donor']) else ''
                }

        cache.add_sample_sequences(sample_id, row_iterator())


def main():
    collapsed_files = snakemake.input.collapsed
    unique_seqs_path = Path(snakemake.output.unique_seqs)
    seq_to_samples_path = Path(snakemake.output.seq_to_samples)
    cache_path = Path(snakemake.output.sqlite_cache)
    log_path = Path(str(snakemake.log))

    log_path.parent.mkdir(parents=True, exist_ok=True)
    cache_path.parent.mkdir(parents=True, exist_ok=True)

    print(f"Creating SQLite cache at: {cache_path}", file=sys.stderr)

    with SequenceCache(cache_path, mode='create') as cache:
        # Process each sample's collapsed sequences
        for i, collapsed_file in enumerate(collapsed_files):
            # Extract sample_id from path
            sample_id = Path(collapsed_file).parent.parent.name

            print(f"Processing sample {i+1}/{len(collapsed_files)}: {sample_id}",
                  file=sys.stderr)

            process_collapsed_file(cache, Path(collapsed_file), sample_id)

        # Get statistics
        stats = cache.get_statistics()

        print(f"\nExporting to TSV files (fast bulk export)...", file=sys.stderr)

        # FAST BULK EXPORT using pandas read_sql (much faster than row-by-row)
        import sqlite3
        conn = sqlite3.connect(str(cache_path))

        # Export unique pairs - direct SQL to DataFrame to gzip
        print(f"  Exporting unique pairs...", file=sys.stderr)
        df_pairs = pd.read_sql_query("SELECT seq_hash, r1_sequence, r2_sequence, total_count FROM unique_pairs", conn)
        df_pairs.to_csv(unique_seqs_path, sep='\t', index=False, compression='gzip')
        print(f"  Exported {len(df_pairs)} unique pairs", file=sys.stderr)

        # Export sample mappings
        print(f"  Exporting sample mappings...", file=sys.stderr)
        df_mappings = pd.read_sql_query("SELECT seq_hash, sample_id, guide, donor, count FROM sample_mappings", conn)
        df_mappings.to_csv(seq_to_samples_path, sep='\t', index=False, compression='gzip')
        print(f"  Exported {len(df_mappings)} sample mappings", file=sys.stderr)

        conn.close()

    # Calculate derived statistics
    total_reads = stats['total_reads']
    single_end_pct = 100 * stats['single_end_reads'] / total_reads if total_reads > 0 else 0
    paired_end_pct = 100 * stats['paired_end_reads'] / total_reads if total_reads > 0 else 0

    n_unique_global = stats['n_unique_pairs']
    n_mappings = stats['n_sample_mappings']
    compression = n_mappings / n_unique_global if n_unique_global > 0 else 1

    log_msg = (
        f"Level 1 Cache - Global Unique R1/R2 Pairs (SQLite):\n"
        f"  Number of samples: {stats['n_samples']}\n"
        f"  Total sample sequence entries: {n_mappings}\n"
        f"  \n"
        f"  Sequencing modes:\n"
        f"    Single-end (R1 only): {stats['single_end_reads']} reads ({single_end_pct:.1f}%)\n"
        f"    Paired-end (R1 + R2): {stats['paired_end_reads']} reads ({paired_end_pct:.1f}%)\n"
        f"  \n"
        f"  Level 1 (alignment cache):\n"
        f"    Unique R1/R2 pairs to align: {n_unique_global}\n"
        f"    Compression vs per-sample: {compression:.1f}x\n"
        f"  \n"
        f"  Level 2 (classification cache):\n"
        f"    Unique (pair, guide, donor) combos: {stats['n_level2_combos']}\n"
        f"    Additional classifications needed: {stats['n_level2_combos'] - n_unique_global}\n"
        f"  \n"
        f"  SAVINGS:\n"
        f"    Alignments saved: {n_mappings - n_unique_global} ({100*(1 - n_unique_global/max(n_mappings, 1)):.1f}%)\n"
        f"  \n"
        f"  SQLite cache: {cache_path}\n"
        f"  ⚠️  PAIRED-END MODE: R1 and R2 will be aligned independently\n"
    )
    log_path.write_text(log_msg)
    print(log_msg, file=sys.stderr)


if __name__ == '__main__':
    main()
