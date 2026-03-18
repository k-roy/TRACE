#!/usr/bin/env python3
"""
Optimized alignment cache with pre-indexing for fast random access.

Author: Kevin R. Roy
Email: kevinrjroy@gmail.com

Key optimizations:
1. PRE-INDEXING: Build index once (seq_hash -> byte offset), reuse forever
2. RANDOM ACCESS: Jump directly to needed rows instead of scanning 1.8GB
3. PERSISTENT CACHE: Store index on disk, rebuild only if source changes
4. PARALLEL INDEXING: Build index using multiple processes
"""

import gzip
import json
import hashlib
import multiprocessing as mp
import os
import sqlite3
import struct
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple
import mmap


def get_file_hash(filepath: Path) -> str:
    """Get hash of file metadata (size + mtime) for cache invalidation."""
    stat = filepath.stat()
    return hashlib.md5(f"{stat.st_size}:{stat.st_mtime}".encode()).hexdigest()[:16]


@dataclass
class IndexEntry:
    """Index entry for a sequence hash."""
    offset: int  # Byte offset in uncompressed data
    length: int  # Length of the line


class AlignmentCacheIndex:
    """
    Pre-indexed alignment cache for O(1) lookups.

    Instead of scanning 1.8GB to find sequences, we:
    1. Build an index: seq_hash -> (byte_offset, line_length)
    2. Store index in SQLite for fast key-value lookups
    3. Use random access to read only the rows we need
    """

    INDEX_VERSION = 2  # Bump when index format changes

    def __init__(self, alignment_cache_path: Path, cache_dir: Optional[Path] = None):
        """
        Initialize with path to alignment cache TSV.gz.

        Args:
            alignment_cache_path: Path to level1_alignment_cache.tsv.gz
            cache_dir: Directory to store index (default: same as alignment cache)
        """
        self.source_path = Path(alignment_cache_path)
        self.cache_dir = cache_dir or self.source_path.parent

        # Index stored as SQLite for fast key lookups
        file_hash = get_file_hash(self.source_path)
        self.index_path = self.cache_dir / f"alignment_index_v{self.INDEX_VERSION}_{file_hash}.db"

        # Decompressed cache for random access
        self.decompressed_path = self.cache_dir / f"alignment_cache_{file_hash}.tsv"

        self.conn = None
        self.header = None
        self.columns = None

    def ensure_index_exists(self, n_workers: int = 8) -> bool:
        """
        Ensure index exists, building it if necessary.

        Returns True if index was already valid, False if rebuilt.
        """
        if self._validate_index():
            print(f"  Using existing index: {self.index_path.name}")
            self._open_index()
            return True

        print(f"  Building new index (one-time cost)...")
        self._build_index(n_workers)
        self._open_index()
        return False

    def _validate_index(self) -> bool:
        """Check if index exists and is valid."""
        if not self.index_path.exists():
            return False
        if not self.decompressed_path.exists():
            return False

        # Check version and source hash
        try:
            conn = sqlite3.connect(self.index_path)
            cursor = conn.execute("SELECT value FROM metadata WHERE key = 'source_hash'")
            row = cursor.fetchone()
            conn.close()

            if row is None:
                return False

            stored_hash = row[0]
            current_hash = get_file_hash(self.source_path)
            return stored_hash == current_hash
        except:
            return False

    def _open_index(self):
        """Open the index for querying."""
        self.conn = sqlite3.connect(self.index_path)

        # Load header info
        cursor = self.conn.execute("SELECT value FROM metadata WHERE key = 'header'")
        self.header = cursor.fetchone()[0]
        self.columns = self.header.split('\t')

        # Open decompressed file for random access
        self.data_file = open(self.decompressed_path, 'r')

    def _build_index(self, n_workers: int = 8):
        """Build the index from the gzipped alignment cache."""
        print(f"    Step 1/3: Decompressing {self.source_path.name}...")

        # Decompress to flat file for random access
        total_bytes = 0
        with gzip.open(self.source_path, 'rt') as f_in:
            with open(self.decompressed_path, 'w') as f_out:
                for line in f_in:
                    f_out.write(line)
                    total_bytes += len(line.encode('utf-8'))

        print(f"    Decompressed: {total_bytes / 1e9:.2f} GB")

        print(f"    Step 2/3: Building index...")

        # Create SQLite index
        if self.index_path.exists():
            self.index_path.unlink()

        conn = sqlite3.connect(self.index_path)
        conn.execute("PRAGMA journal_mode = WAL")
        conn.execute("PRAGMA synchronous = NORMAL")
        conn.execute("PRAGMA cache_size = -64000")  # 64MB cache

        # Create tables
        conn.execute("""
            CREATE TABLE metadata (
                key TEXT PRIMARY KEY,
                value TEXT
            )
        """)
        conn.execute("""
            CREATE TABLE seq_index (
                seq_hash TEXT PRIMARY KEY,
                byte_offset INTEGER,
                line_length INTEGER
            )
        """)

        # Read and index
        with open(self.decompressed_path, 'r') as f:
            # Read header
            header = f.readline().strip()
            header_offset = len(header.encode('utf-8')) + 1  # +1 for newline

            conn.execute("INSERT INTO metadata VALUES (?, ?)", ('header', header))
            conn.execute("INSERT INTO metadata VALUES (?, ?)", ('source_hash', get_file_hash(self.source_path)))

            # Index all rows
            current_offset = header_offset
            batch = []
            indexed = 0

            for line in f:
                line_bytes = len(line.encode('utf-8'))
                seq_hash = line.split('\t', 1)[0]

                batch.append((seq_hash, current_offset, line_bytes))
                current_offset += line_bytes

                if len(batch) >= 50000:
                    conn.executemany(
                        "INSERT OR REPLACE INTO seq_index VALUES (?, ?, ?)",
                        batch
                    )
                    indexed += len(batch)
                    if indexed % 500000 == 0:
                        print(f"      Indexed {indexed:,} sequences...")
                    batch = []

            # Final batch
            if batch:
                conn.executemany(
                    "INSERT OR REPLACE INTO seq_index VALUES (?, ?, ?)",
                    batch
                )
                indexed += len(batch)

        print(f"    Step 3/3: Creating index on seq_hash...")
        # Index is already on primary key, just commit
        conn.commit()
        conn.close()

        print(f"    Index complete: {indexed:,} sequences indexed")

    def get(self, seq_hash: str) -> Optional[Dict]:
        """
        Get alignment data for a sequence hash.

        O(1) lookup using pre-built index.
        """
        cursor = self.conn.execute(
            "SELECT byte_offset, line_length FROM seq_index WHERE seq_hash = ?",
            (seq_hash,)
        )
        row = cursor.fetchone()

        if row is None:
            return None

        offset, length = row

        # Random access read
        self.data_file.seek(offset)
        line = self.data_file.read(length).strip()

        # Parse the line
        values = line.split('\t')
        result = {}

        for i, col in enumerate(self.columns):
            if i >= len(values):
                result[col] = ''
                continue

            val = values[i]

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

    def get_batch(self, seq_hashes: List[str]) -> Dict[str, Dict]:
        """
        Get alignment data for multiple sequence hashes efficiently.

        Groups lookups by offset to minimize seeks.
        """
        if not seq_hashes:
            return {}

        # Get all offsets in one query
        placeholders = ','.join(['?'] * len(seq_hashes))
        cursor = self.conn.execute(
            f"SELECT seq_hash, byte_offset, line_length FROM seq_index WHERE seq_hash IN ({placeholders})",
            seq_hashes
        )

        # Sort by offset for sequential reading
        entries = sorted(cursor.fetchall(), key=lambda x: x[1])

        results = {}
        for seq_hash, offset, length in entries:
            self.data_file.seek(offset)
            line = self.data_file.read(length).strip()

            values = line.split('\t')
            result = {}

            for i, col in enumerate(self.columns):
                if i >= len(values):
                    result[col] = ''
                    continue

                val = values[i]

                if col in ('r1_ref_start', 'r1_ref_end', 'r2_ref_start', 'r2_ref_end'):
                    try:
                        result[col] = int(val) if val and val != 'nan' else -1
                    except (ValueError, TypeError):
                        result[col] = -1
                elif col in ('r1_is_mapped', 'r2_is_mapped', 'is_proper_pair'):
                    result[col] = val == 'True' or val == '1' or val == True
                else:
                    result[col] = val if val != 'nan' else ''

            results[seq_hash] = result

        return results

    def __contains__(self, seq_hash: str) -> bool:
        """Check if sequence hash exists in index."""
        cursor = self.conn.execute(
            "SELECT 1 FROM seq_index WHERE seq_hash = ? LIMIT 1",
            (seq_hash,)
        )
        return cursor.fetchone() is not None

    def count(self) -> int:
        """Get total number of indexed sequences."""
        cursor = self.conn.execute("SELECT COUNT(*) FROM seq_index")
        return cursor.fetchone()[0]

    def close(self):
        """Close index and data file."""
        if self.conn:
            self.conn.close()
        if hasattr(self, 'data_file') and self.data_file:
            self.data_file.close()


class FastAlignmentLookup:
    """
    Drop-in replacement for SQLiteAlignmentLookup with pre-indexing.

    Key difference: Uses pre-built index for O(1) lookups instead of
    loading relevant hashes into a temp SQLite database every run.
    """

    def __init__(self, alignment_cache_path: Path, relevant_hashes: Set[str],
                 cache_dir: Optional[Path] = None):
        """
        Initialize fast alignment lookup.

        Args:
            alignment_cache_path: Path to level1_alignment_cache.tsv.gz
            relevant_hashes: Set of seq_hashes we need (used for counting, not filtering)
            cache_dir: Directory for index cache (default: same as alignment cache)
        """
        self.index = AlignmentCacheIndex(alignment_cache_path, cache_dir)
        was_cached = self.index.ensure_index_exists()

        self.relevant_hashes = relevant_hashes
        self.loaded_count = len(relevant_hashes)  # For compatibility

        if was_cached:
            print(f"  Index ready (cached): {self.index.count():,} total sequences")
        else:
            print(f"  Index ready (built): {self.index.count():,} total sequences")

    def __contains__(self, seq_hash: str) -> bool:
        return seq_hash in self.index

    def __getitem__(self, seq_hash: str) -> Dict:
        result = self.index.get(seq_hash)
        if result is None:
            raise KeyError(seq_hash)
        return result

    def get(self, seq_hash: str, default=None) -> Optional[Dict]:
        return self.index.get(seq_hash) or default

    def get_batch(self, seq_hashes: List[str]) -> Dict[str, Dict]:
        return self.index.get_batch(seq_hashes)

    def close(self):
        self.index.close()


if __name__ == '__main__':
    import argparse
    import time

    parser = argparse.ArgumentParser(description='Build or test alignment cache index')
    parser.add_argument('alignment_cache', type=Path, help='Path to alignment cache TSV.gz')
    parser.add_argument('--cache-dir', type=Path, help='Directory for index cache')
    parser.add_argument('--test', action='store_true', help='Run lookup tests')
    parser.add_argument('--benchmark', type=int, help='Benchmark N random lookups')

    args = parser.parse_args()

    print(f"Alignment cache: {args.alignment_cache}")

    # Build/load index
    start = time.time()
    index = AlignmentCacheIndex(args.alignment_cache, args.cache_dir)
    index.ensure_index_exists()
    print(f"Index ready in {time.time() - start:.1f}s")
    print(f"Total sequences: {index.count():,}")

    if args.test:
        # Test with a few lookups
        print("\nTesting lookups...")
        cursor = index.conn.execute("SELECT seq_hash FROM seq_index LIMIT 5")
        test_hashes = [row[0] for row in cursor.fetchall()]

        for h in test_hashes:
            start = time.time()
            result = index.get(h)
            elapsed = (time.time() - start) * 1000
            print(f"  {h[:20]}... -> {len(result)} fields ({elapsed:.2f}ms)")

    if args.benchmark:
        import random

        # Get random sample of hashes
        cursor = index.conn.execute(f"SELECT seq_hash FROM seq_index ORDER BY RANDOM() LIMIT {args.benchmark}")
        test_hashes = [row[0] for row in cursor.fetchall()]

        print(f"\nBenchmarking {args.benchmark} random lookups...")

        # Individual lookups
        start = time.time()
        for h in test_hashes:
            _ = index.get(h)
        individual_time = time.time() - start
        print(f"  Individual: {individual_time:.2f}s ({args.benchmark / individual_time:.0f} lookups/sec)")

        # Batch lookup
        start = time.time()
        _ = index.get_batch(test_hashes)
        batch_time = time.time() - start
        print(f"  Batch: {batch_time:.2f}s ({args.benchmark / batch_time:.0f} lookups/sec)")

    index.close()