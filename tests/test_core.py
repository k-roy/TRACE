"""Tests for trace.core modules."""

import pytest
from trace_crispr.core.scoring import (
    AlignmentScore,
    count_mismatches_from_md,
    DeduplicationSignature,
)
from trace_crispr.core.classification import (
    EditingOutcome,
    ClassificationResult,
    get_hdr_signature_positions,
)
from trace_crispr.core.kmer import KmerClassifier


class TestAlignmentScore:
    """Test AlignmentScore class."""

    def test_penalty_score_perfect_match(self):
        """Test penalty score for perfect match."""
        score = AlignmentScore(
            aligner="bwa",
            mismatches=0,
            soft_clips=0,
            is_aligned=True,
        )
        assert score.penalty_score == 0

    def test_penalty_score_with_mismatches(self):
        """Test penalty score with mismatches."""
        score = AlignmentScore(
            aligner="bwa",
            mismatches=5,
            soft_clips=0,
            is_aligned=True,
        )
        assert score.penalty_score == 5

    def test_penalty_score_with_soft_clips(self):
        """Test penalty score with soft clips."""
        score = AlignmentScore(
            aligner="bwa",
            mismatches=0,
            soft_clips=20,
            is_aligned=True,
        )
        assert score.penalty_score == 20

    def test_penalty_score_combined(self):
        """Test penalty score with mismatches and soft clips."""
        score = AlignmentScore(
            aligner="bwa",
            mismatches=3,
            soft_clips=10,
            is_aligned=True,
        )
        assert score.penalty_score == 13  # 3 + 10

    def test_penalty_score_unaligned(self):
        """Test penalty score for unaligned reads is infinity."""
        score = AlignmentScore(aligner="bwa", is_aligned=False)
        assert score.penalty_score == float('inf')

    def test_priority_bwa(self):
        """Test BWA has highest priority."""
        score = AlignmentScore(aligner="bwa")
        assert score.priority == 0

    def test_priority_bbmap(self):
        """Test BBMap has second priority."""
        score = AlignmentScore(aligner="bbmap")
        assert score.priority == 1

    def test_priority_minimap2(self):
        """Test minimap2 has third priority."""
        score = AlignmentScore(aligner="minimap2")
        assert score.priority == 2

    def test_comparison_lower_penalty_wins(self):
        """Test alignment with lower penalty wins."""
        score1 = AlignmentScore(aligner="minimap2", mismatches=1, soft_clips=0, is_aligned=True)
        score2 = AlignmentScore(aligner="bwa", mismatches=5, soft_clips=0, is_aligned=True)
        assert score1 < score2

    def test_comparison_tiebreaker_priority(self):
        """Test tie-breaker uses priority."""
        score1 = AlignmentScore(aligner="bwa", mismatches=2, soft_clips=0, is_aligned=True)
        score2 = AlignmentScore(aligner="minimap2", mismatches=2, soft_clips=0, is_aligned=True)
        assert score1 < score2  # BWA wins tie-breaker


class TestMDTagParsing:
    """Test MD tag parsing."""

    def test_count_mismatches_no_mismatches(self):
        """Test MD tag with no mismatches."""
        assert count_mismatches_from_md("100") == 0

    def test_count_mismatches_single(self):
        """Test MD tag with single mismatch."""
        assert count_mismatches_from_md("50A50") == 1

    def test_count_mismatches_multiple(self):
        """Test MD tag with multiple mismatches."""
        assert count_mismatches_from_md("30A20T50") == 2

    def test_count_mismatches_with_deletion(self):
        """Test MD tag with deletion (should not count as mismatch)."""
        assert count_mismatches_from_md("50^ATG50") == 0

    def test_count_mismatches_complex(self):
        """Test complex MD tag."""
        # 30 matches, mismatch A, 10 matches, deletion ATG, 20 matches, mismatch C
        assert count_mismatches_from_md("30A10^ATG20C40") == 2


class TestDeduplicationSignature:
    """Test deduplication signature."""

    def test_signature_equality(self):
        """Test that identical signatures are equal."""
        sig1 = DeduplicationSignature(
            r1_ref_start=0,
            r1_ref_end=100,
            r2_ref_start=50,
            r2_ref_end=150,
            r1_cigar="100M",
            r2_cigar="100M",
        )
        sig2 = DeduplicationSignature(
            r1_ref_start=0,
            r1_ref_end=100,
            r2_ref_start=50,
            r2_ref_end=150,
            r1_cigar="100M",
            r2_cigar="100M",
        )
        assert sig1 == sig2

    def test_signature_hash(self):
        """Test that identical signatures have same hash."""
        sig1 = DeduplicationSignature(0, 100, 50, 150, "100M", "100M")
        sig2 = DeduplicationSignature(0, 100, 50, 150, "100M", "100M")
        assert hash(sig1) == hash(sig2)

    def test_signature_in_set(self):
        """Test that signatures work in sets for deduplication."""
        sig1 = DeduplicationSignature(0, 100, 50, 150, "100M", "100M")
        sig2 = DeduplicationSignature(0, 100, 50, 150, "100M", "100M")  # Duplicate
        sig3 = DeduplicationSignature(1, 101, 51, 151, "100M", "100M")  # Different

        seen = {sig1}
        assert sig2 in seen
        assert sig3 not in seen


class TestEditingOutcome:
    """Test EditingOutcome enum."""

    def test_outcome_values(self):
        """Test editing outcome enum values."""
        assert EditingOutcome.WILD_TYPE.value == "wild_type"
        assert EditingOutcome.HDR_PERFECT.value == "hdr_perfect"
        assert EditingOutcome.HDR_IMPERFECT.value == "hdr_imperfect"
        assert EditingOutcome.NHEJ_INSERTION.value == "nhej_insertion"
        assert EditingOutcome.NHEJ_DELETION.value == "nhej_deletion"
        assert EditingOutcome.LARGE_DELETION.value == "large_deletion"
        assert EditingOutcome.UNCLASSIFIED.value == "unclassified"


class TestHDRSignature:
    """Test HDR signature detection."""

    def test_get_hdr_signature_positions(self):
        """Test detecting HDR signature positions."""
        reference = "AAACCCGGG"
        hdr_template = "AAATTTGGG"

        positions = get_hdr_signature_positions(reference, hdr_template)

        assert len(positions) == 3
        # Positions 3, 4, 5 differ (CCC -> TTT)
        assert (3, 'C', 'T') in positions
        assert (4, 'C', 'T') in positions
        assert (5, 'C', 'T') in positions

    def test_no_differences(self):
        """Test identical sequences have no signature."""
        sequence = "AAACCCGGG"
        positions = get_hdr_signature_positions(sequence, sequence)
        assert len(positions) == 0


class TestKmerClassifier:
    """Test k-mer classifier."""

    def test_classifier_initialization(self):
        """Test k-mer classifier initialization."""
        classifier = KmerClassifier.from_sequences(
            reference="AAACCCGGGTTTATATATAT",
            hdr_template="AAATTTGGGTTTATATATAT",
            edit_positions=[3, 4, 5],
            kmer_size=5,
        )

        assert classifier.kmer_size == 5
        assert len(classifier.wt_kmers) > 0
        assert len(classifier.hdr_kmers) > 0

    def test_classify_wt_read(self):
        """Test classifying a WT read."""
        classifier = KmerClassifier.from_sequences(
            reference="AAACCCGGGTTTATATATAT",
            hdr_template="AAATTTGGGTTTATATATAT",
            edit_positions=[3, 4, 5],
            kmer_size=5,
        )

        # Read matching WT (contains CCC)
        result = classifier.classify_sequence("AAACCCGGGTT")
        # Result is an enum; check its value
        assert result.value in ("wt", "ambiguous", "unknown")

    def test_classify_hdr_read(self):
        """Test classifying an HDR read."""
        classifier = KmerClassifier.from_sequences(
            reference="AAACCCGGGTTTATATATAT",
            hdr_template="AAATTTGGGTTTATATATAT",
            edit_positions=[3, 4, 5],
            kmer_size=5,
        )

        # Read matching HDR (contains TTT)
        result = classifier.classify_sequence("AAATTTGGGTT")
        # Result is an enum; check its value
        assert result.value in ("hdr", "ambiguous", "unknown")


class TestClassificationResult:
    """Test ClassificationResult dataclass."""

    def test_classification_result_creation(self):
        """Test creating a classification result."""
        result = ClassificationResult(
            outcome=EditingOutcome.HDR_PERFECT,
            confidence=0.95,
            details={"num_hdr_positions": 2},
        )

        assert result.outcome == EditingOutcome.HDR_PERFECT
        assert result.confidence == 0.95
        assert result.details["num_hdr_positions"] == 2


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
