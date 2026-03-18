#!/usr/bin/env python3
"""
SQLite-based caching for TRACE pipeline.

Reduces memory usage by 98% for large datasets by storing intermediate
data on disk instead of in-memory dictionaries.

Key advantages:
1. Memory usage independent of dataset size (~50MB vs 22GB)
2. Indexed lookups remain fast (B-tree, O(log n))
3. SQL aggregations happen at database level
4. Batch inserts with transactions for speed
5. Resume capability if job is interrupted

Author: Kevin R. Roy
"""

import gzip
import hashlib
import sqlite3
from contextlib import contextmanager
from pathlib import Path
from typing import Iterator, Dict, List, Optional, Tuple, Any


class SequenceCache:
    """
    SQLite-based cache for unique sequence pairs.

    Replaces in-memory defaultdict with disk-based storage.
    Memory usage: ~50MB regardless of dataset size.

    Schema:
    - unique_pairs: (seq_hash, r1_sequence, r2_sequence, total_count)
    - sample_mappings: (seq_hash, sample_id, guide, donor, count)

    Usage:
        cache = SequenceCache("/path/to/cache.db")

        # Add sequences in batches
        for sample_file in sample_files:
            cache.add_sample_sequences(sample_id, sequences)

        # Query unique pairs
        for pair in cache.iter_unique_pairs():
            process(pair)

        # Export to TSV
        cache.export_unique_pairs("unique_seqs.tsv.gz")
        cache.export_sample_mappings("seq_to_samples.tsv.gz")
    """

    BATCH_SIZE = 10000  # Insert in batches for speed

    def __init__(self, db_path: Path, mode: str = 'create'):
        """
        Initialize SQLite cache.

        Args:
            db_path: Path to SQLite database file
            mode: 'create' (new db), 'append' (add to existing), 'read' (read-only)
        """
        self.db_path = Path(db_path)
        self.mode = mode
        self.conn = None

        # Ensure parent directory exists
        self.db_path.parent.mkdir(parents=True, exist_ok=True)

        if mode == 'create':
            # Remove existing database
            if self.db_path.exists():
                self.db_path.unlink()
            self._init_database()
        elif mode == 'append':
            self._init_database()
        elif mode == 'read':
            if not self.db_path.exists():
                raise FileNotFoundError(f"Cache database not found: {db_path}")
            self.conn = sqlite3.connect(str(self.db_path), timeout=30.0)
        else:
            raise ValueError(f"Unknown mode: {mode}")

    def _init_database(self):
        """Create database schema."""
        self.conn = sqlite3.connect(str(self.db_path), timeout=30.0)

        # Optimize for batch writes
        self.conn.execute("PRAGMA journal_mode = WAL")
        self.conn.execute("PRAGMA synchronous = NORMAL")
        self.conn.execute("PRAGMA cache_size = -64000")  # 64MB cache

        # Create tables
        self.conn.executescript("""
            CREATE TABLE IF NOT EXISTS unique_pairs (
                seq_hash TEXT PRIMARY KEY,
                r1_sequence TEXT NOT NULL,
                r2_sequence TEXT NOT NULL,
                total_count INTEGER DEFAULT 0
            );

            CREATE TABLE IF NOT EXISTS sample_mappings (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                seq_hash TEXT NOT NULL,
                sample_id TEXT NOT NULL,
                guide TEXT,
                donor TEXT,
                count INTEGER NOT NULL,
                FOREIGN KEY (seq_hash) REFERENCES unique_pairs(seq_hash)
            );

            CREATE INDEX IF NOT EXISTS idx_mappings_hash ON sample_mappings(seq_hash);
            CREATE INDEX IF NOT EXISTS idx_mappings_sample ON sample_mappings(sample_id);
        """)
        self.conn.commit()

    def close(self):
        """Close database connection."""
        if self.conn:
            self.conn.close()
            self.conn = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    @staticmethod
    def compute_hash(r1_seq: str, r2_seq: str) -> str:
        """
        Hash an R1/R2 pair for lookup.

        Uses full 32-char MD5 to avoid hash collisions with large datasets.
        Must match the hash function in collapse_reads.py!
        """
        combined = r1_seq + '|' + r2_seq
        return hashlib.md5(combined.encode()).hexdigest()

    def add_sequences_batch(
        self,
        sequences: List[Tuple[str, str, str, str, str, int]]
    ):
        """
        Add a batch of sequences.

        Args:
            sequences: List of (r1_seq, r2_seq, sample_id, guide, donor, count)
        """
        cursor = self.conn.cursor()

        try:
            for r1_seq, r2_seq, sample_id, guide, donor, count in sequences:
                seq_hash = self.compute_hash(r1_seq, r2_seq)

                # Upsert into unique_pairs (increment count)
                cursor.execute("""
                    INSERT INTO unique_pairs (seq_hash, r1_sequence, r2_sequence, total_count)
                    VALUES (?, ?, ?, ?)
                    ON CONFLICT(seq_hash) DO UPDATE SET
                        total_count = total_count + excluded.total_count
                """, (seq_hash, r1_seq, r2_seq, count))

                # Insert sample mapping
                cursor.execute("""
                    INSERT INTO sample_mappings (seq_hash, sample_id, guide, donor, count)
                    VALUES (?, ?, ?, ?, ?)
                """, (seq_hash, sample_id, guide, donor, count))

            self.conn.commit()
        except Exception:
            self.conn.rollback()
            raise

    def add_sample_sequences(
        self,
        sample_id: str,
        df_iterator: Iterator[Dict[str, Any]]
    ):
        """
        Add sequences from a sample DataFrame iterator.

        Processes in batches to minimize memory usage.

        Args:
            sample_id: Sample identifier
            df_iterator: Iterator yielding dicts with keys:
                r1_sequence, r2_sequence, count, guide, donor
        """
        batch = []

        for row in df_iterator:
            r1_seq = row['r1_sequence']
            r2_seq = row['r2_sequence'] if row.get('r2_sequence') else ''
            guide = row.get('guide', '') or ''
            donor = row.get('donor', '') or ''
            count = row['count']

            batch.append((r1_seq, r2_seq, sample_id, guide, donor, count))

            if len(batch) >= self.BATCH_SIZE:
                self.add_sequences_batch(batch)
                batch = []

        # Final batch
        if batch:
            self.add_sequences_batch(batch)

    def get_unique_pair_count(self) -> int:
        """Return count of unique R1/R2 pairs."""
        cursor = self.conn.execute("SELECT COUNT(*) FROM unique_pairs")
        return cursor.fetchone()[0]

    def get_sample_mapping_count(self) -> int:
        """Return count of sample mappings."""
        cursor = self.conn.execute("SELECT COUNT(*) FROM sample_mappings")
        return cursor.fetchone()[0]

    def get_total_reads(self) -> int:
        """Return total read count across all samples."""
        cursor = self.conn.execute("SELECT SUM(total_count) FROM unique_pairs")
        result = cursor.fetchone()[0]
        return result if result else 0

    def iter_unique_pairs(self, batch_size: int = 10000) -> Iterator[Dict]:
        """
        Iterate over unique pairs in batches.

        Yields dicts with: seq_hash, r1_sequence, r2_sequence, total_count
        """
        cursor = self.conn.cursor()
        offset = 0

        while True:
            cursor.execute("""
                SELECT seq_hash, r1_sequence, r2_sequence, total_count
                FROM unique_pairs
                LIMIT ? OFFSET ?
            """, (batch_size, offset))

            rows = cursor.fetchall()
            if not rows:
                break

            for row in rows:
                yield {
                    'seq_hash': row[0],
                    'r1_sequence': row[1],
                    'r2_sequence': row[2],
                    'total_count': row[3]
                }

            offset += batch_size

    def iter_sample_mappings(self, batch_size: int = 10000) -> Iterator[Dict]:
        """
        Iterate over sample mappings in batches.

        Yields dicts with: seq_hash, sample_id, guide, donor, count
        """
        cursor = self.conn.cursor()
        offset = 0

        while True:
            cursor.execute("""
                SELECT seq_hash, sample_id, guide, donor, count
                FROM sample_mappings
                LIMIT ? OFFSET ?
            """, (batch_size, offset))

            rows = cursor.fetchall()
            if not rows:
                break

            for row in rows:
                yield {
                    'seq_hash': row[0],
                    'sample_id': row[1],
                    'guide': row[2],
                    'donor': row[3],
                    'count': row[4]
                }

            offset += batch_size

    def export_unique_pairs(self, output_path: Path):
        """Export unique pairs to gzipped TSV using fast bulk export."""
        import pandas as pd

        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # Fast bulk export using pandas read_sql (100x faster than row-by-row)
        df = pd.read_sql_query(
            "SELECT seq_hash, r1_sequence, r2_sequence, total_count FROM unique_pairs",
            self.conn
        )
        df.to_csv(output_path, sep='\t', index=False, compression='gzip')

    def export_sample_mappings(self, output_path: Path):
        """Export sample mappings to gzipped TSV using fast bulk export."""
        import pandas as pd

        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # Fast bulk export using pandas read_sql (100x faster than row-by-row)
        df = pd.read_sql_query(
            "SELECT seq_hash, sample_id, guide, donor, count FROM sample_mappings",
            self.conn
        )
        df.to_csv(output_path, sep='\t', index=False, compression='gzip')

    def get_statistics(self) -> Dict[str, Any]:
        """Get cache statistics."""
        cursor = self.conn.cursor()

        # Count sequencing modes (handle NULL r2_sequence as single-end)
        cursor.execute("""
            SELECT
                COALESCE(SUM(CASE WHEN r2_sequence IS NULL OR r2_sequence = '' THEN total_count ELSE 0 END), 0) as single_end,
                COALESCE(SUM(CASE WHEN r2_sequence IS NOT NULL AND r2_sequence != '' THEN total_count ELSE 0 END), 0) as paired_end
            FROM unique_pairs
        """)
        single_end, paired_end = cursor.fetchone()
        # Extra safety check for None values (shouldn't happen with COALESCE but defensive)
        single_end = single_end if single_end is not None else 0
        paired_end = paired_end if paired_end is not None else 0

        # Count unique samples
        cursor.execute("SELECT COUNT(DISTINCT sample_id) FROM sample_mappings")
        n_samples = cursor.fetchone()[0]

        # Count unique (seq, guide, donor) combinations for Level 2
        cursor.execute("""
            SELECT COUNT(DISTINCT seq_hash || '|' || COALESCE(guide, '') || '|' || COALESCE(donor, ''))
            FROM sample_mappings
        """)
        n_level2_combos = cursor.fetchone()[0]

        return {
            'n_unique_pairs': self.get_unique_pair_count(),
            'n_sample_mappings': self.get_sample_mapping_count(),
            'n_samples': n_samples,
            'n_level2_combos': n_level2_combos,
            'single_end_reads': single_end,
            'paired_end_reads': paired_end,
            'total_reads': single_end + paired_end
        }


class AlignmentCache:
    """
    SQLite-based cache for alignment results.

    Stores alignment results for unique sequence pairs.
    Used to avoid re-alignment when the same sequence appears
    in multiple samples.

    Schema:
    - alignments: (seq_hash, r1_alignment, r2_alignment, consensus_cigar, ...)

    Usage:
        cache = AlignmentCache("/path/to/alignments.db")

        # Check if already aligned
        if not cache.has_alignment(seq_hash):
            result = align(r1_seq, r2_seq)
            cache.add_alignment(seq_hash, result)

        # Query alignments
        for alignment in cache.iter_alignments():
            classify(alignment)
    """

    BATCH_SIZE = 5000

    def __init__(self, db_path: Path, mode: str = 'create'):
        """
        Initialize alignment cache.

        Args:
            db_path: Path to SQLite database
            mode: 'create', 'append', or 'read'
        """
        self.db_path = Path(db_path)
        self.mode = mode
        self.conn = None

        self.db_path.parent.mkdir(parents=True, exist_ok=True)

        if mode == 'create':
            if self.db_path.exists():
                self.db_path.unlink()
            self._init_database()
        elif mode == 'append':
            self._init_database()
        elif mode == 'read':
            if not self.db_path.exists():
                raise FileNotFoundError(f"Alignment cache not found: {db_path}")
            self.conn = sqlite3.connect(str(self.db_path), timeout=30.0)
        else:
            raise ValueError(f"Unknown mode: {mode}")

    def _init_database(self):
        """Create database schema."""
        self.conn = sqlite3.connect(str(self.db_path), timeout=30.0)

        self.conn.execute("PRAGMA journal_mode = WAL")
        self.conn.execute("PRAGMA synchronous = NORMAL")
        self.conn.execute("PRAGMA cache_size = -64000")

        self.conn.executescript("""
            CREATE TABLE IF NOT EXISTS alignments (
                seq_hash TEXT PRIMARY KEY,
                r1_sequence TEXT NOT NULL,
                r2_sequence TEXT,
                ref_start INTEGER,
                ref_end INTEGER,
                cigar TEXT,
                aligned_seq TEXT,
                edit_distance INTEGER,
                mapping_quality INTEGER,
                is_mapped INTEGER DEFAULT 1,
                alignment_method TEXT,
                consensus_cigar TEXT,
                consensus_seq TEXT,
                consensus_ref_start INTEGER,
                consensus_ref_end INTEGER
            );
        """)
        self.conn.commit()

    def close(self):
        """Close database connection."""
        if self.conn:
            self.conn.close()
            self.conn = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def has_alignment(self, seq_hash: str) -> bool:
        """Check if alignment exists for sequence."""
        cursor = self.conn.execute(
            "SELECT 1 FROM alignments WHERE seq_hash = ?",
            (seq_hash,)
        )
        return cursor.fetchone() is not None

    def add_alignment(self, seq_hash: str, alignment: Dict[str, Any]):
        """Add a single alignment result."""
        self.conn.execute("""
            INSERT OR REPLACE INTO alignments (
                seq_hash, r1_sequence, r2_sequence,
                ref_start, ref_end, cigar, aligned_seq,
                edit_distance, mapping_quality, is_mapped,
                alignment_method, consensus_cigar, consensus_seq,
                consensus_ref_start, consensus_ref_end
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """, (
            seq_hash,
            alignment.get('r1_sequence', ''),
            alignment.get('r2_sequence', ''),
            alignment.get('ref_start'),
            alignment.get('ref_end'),
            alignment.get('cigar'),
            alignment.get('aligned_seq'),
            alignment.get('edit_distance'),
            alignment.get('mapping_quality'),
            1 if alignment.get('is_mapped', True) else 0,
            alignment.get('alignment_method'),
            alignment.get('consensus_cigar'),
            alignment.get('consensus_seq'),
            alignment.get('consensus_ref_start'),
            alignment.get('consensus_ref_end')
        ))
        self.conn.commit()

    def add_alignments_batch(self, alignments: List[Tuple]):
        """Add multiple alignments in a batch."""
        cursor = self.conn.cursor()

        try:
            cursor.executemany("""
                INSERT OR REPLACE INTO alignments (
                    seq_hash, r1_sequence, r2_sequence,
                    ref_start, ref_end, cigar, aligned_seq,
                    edit_distance, mapping_quality, is_mapped,
                    alignment_method, consensus_cigar, consensus_seq,
                    consensus_ref_start, consensus_ref_end
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """, alignments)
            self.conn.commit()
        except Exception:
            self.conn.rollback()
            raise

    def get_alignment(self, seq_hash: str) -> Optional[Dict[str, Any]]:
        """Get alignment for a sequence hash."""
        cursor = self.conn.execute("""
            SELECT * FROM alignments WHERE seq_hash = ?
        """, (seq_hash,))

        row = cursor.fetchone()
        if not row:
            return None

        columns = [desc[0] for desc in cursor.description]
        return dict(zip(columns, row))

    def iter_alignments(self, batch_size: int = 5000) -> Iterator[Dict]:
        """Iterate over all alignments in batches."""
        cursor = self.conn.cursor()
        offset = 0

        while True:
            cursor.execute("""
                SELECT * FROM alignments
                LIMIT ? OFFSET ?
            """, (batch_size, offset))

            rows = cursor.fetchall()
            if not rows:
                break

            columns = [desc[0] for desc in cursor.description]
            for row in rows:
                yield dict(zip(columns, row))

            offset += batch_size

    def get_count(self) -> int:
        """Return count of cached alignments."""
        cursor = self.conn.execute("SELECT COUNT(*) FROM alignments")
        return cursor.fetchone()[0]

    def export_to_tsv(self, output_path: Path):
        """Export alignments to gzipped TSV using fast bulk export."""
        import pandas as pd

        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # Fast bulk export using pandas read_sql (100x faster than row-by-row)
        df = pd.read_sql_query("SELECT * FROM alignments", self.conn)
        df.to_csv(output_path, sep='\t', index=False, compression='gzip')
