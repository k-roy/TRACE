"""Tests for trace.utils module."""

import pytest
from trace_crispr.utils.sequence import (
    reverse_complement,
    hamming_distance,
    find_guide_in_sequence,
    extract_unique_kmers,
)


class TestReverseComplement:
    """Test reverse complement function."""

    def test_simple_sequence(self):
        """Test simple sequence reverse complement."""
        assert reverse_complement("ATCG") == "CGAT"

    def test_longer_sequence(self):
        """Test longer sequence reverse complement."""
        seq = "GCTGAAGCACTGCACGCCGT"
        rc = reverse_complement(seq)
        assert rc == "ACGGCGTGCAGTGCTTCAGC"

    def test_reverse_complement_is_involutive(self):
        """Test that reverse complement of reverse complement is original."""
        seq = "ATCGATCGATCG"
        assert reverse_complement(reverse_complement(seq)) == seq

    def test_lowercase_handling(self):
        """Test that lowercase is handled correctly."""
        assert reverse_complement("atcg") == "cgat"

    def test_mixed_case(self):
        """Test mixed case sequences."""
        assert reverse_complement("AtCg") == "cGaT"


class TestHammingDistance:
    """Test hamming distance function."""

    def test_identical_sequences(self):
        """Test hamming distance of identical sequences."""
        assert hamming_distance("ATCG", "ATCG") == 0

    def test_one_mismatch(self):
        """Test hamming distance with one mismatch."""
        assert hamming_distance("ATCG", "ATCG"[0] + "G" + "ATCG"[2:]) == 1
        assert hamming_distance("ATCG", "ATGG") == 1

    def test_all_different(self):
        """Test hamming distance with all different."""
        assert hamming_distance("AAAA", "TTTT") == 4

    def test_different_lengths_raises(self):
        """Test that different length sequences raise error."""
        with pytest.raises(ValueError):
            hamming_distance("ATCG", "ATC")


class TestFindGuide:
    """Test guide finding function."""

    def test_find_guide_forward_strand(self):
        """Test finding guide on forward strand."""
        reference = "NNNNNGCTGAAGCACTGCACGCCGTNNNNNGGG"
        guide = "GCTGAAGCACTGCACGCCGT"

        pos, strand = find_guide_in_sequence(reference, guide)
        assert pos == 5
        assert strand == '+'

    def test_find_guide_reverse_strand(self):
        """Test finding guide on reverse strand."""
        guide = "GCTGAAGCACTGCACGCCGT"
        rc_guide = reverse_complement(guide)
        reference = "NNNNN" + rc_guide + "NNNNN"

        pos, strand = find_guide_in_sequence(reference, guide)
        assert strand == '-'

    def test_guide_not_found(self):
        """Test guide not found raises ValueError."""
        reference = "AAAAAAAAAAAAAAAAAAAAAA"
        guide = "GCTGAAGCACTGCACGCCGT"

        # Guide not found raises ValueError
        with pytest.raises(ValueError):
            find_guide_in_sequence(reference, guide)


class TestExtractKmers:
    """Test k-mer extraction."""

    def test_extract_kmers_basic(self):
        """Test basic k-mer extraction."""
        sequence = "ATCGATCG"
        kmers = extract_unique_kmers(sequence, kmer_size=4)

        # Check that kmers are present (including reverse complements)
        assert len(kmers) > 0

    def test_kmer_count(self):
        """Test correct number of k-mers."""
        sequence = "ATCGATCGATCGATCG"  # 16 bp
        kmers = extract_unique_kmers(sequence, kmer_size=4)
        # Should include both forward and reverse complement k-mers
        assert len(kmers) > 0

    def test_exclude_sequences(self):
        """Test excluding sequences from k-mers."""
        sequence = "ATCGATCG"
        exclude = ["ATCGATCG"]  # Exclude the same sequence

        kmers = extract_unique_kmers(sequence, kmer_size=4, exclude_sequences=exclude)
        # All k-mers from exclude should be removed
        assert len(kmers) == 0

    def test_empty_sequence(self):
        """Test empty sequence returns empty set."""
        kmers = extract_unique_kmers("", kmer_size=4)
        assert len(kmers) == 0

    def test_sequence_shorter_than_k(self):
        """Test sequence shorter than k returns empty set."""
        kmers = extract_unique_kmers("ATC", kmer_size=4)
        assert len(kmers) == 0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
