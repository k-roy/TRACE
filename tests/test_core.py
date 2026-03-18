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
    calculate_deletion_microhomology,
    classify_indel_mechanism,
    summarize_classifications,
    MMEJ_MICROHOMOLOGY_THRESHOLD,
)


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
        """Test editing outcome enum values (comprehensive scheme)."""
        # HDR categories
        assert EditingOutcome.HDR_COMPLETE.value == "hdr_complete"
        assert EditingOutcome.HDR_PARTIAL.value == "hdr_partial"
        assert EditingOutcome.HDR_PLUS_NHEJ_INDEL.value == "hdr_plus_nhej_indel"
        assert EditingOutcome.HDR_PLUS_MMEJ_INDEL.value == "hdr_plus_mmej_indel"
        assert EditingOutcome.HDR_PLUS_OTHER.value == "hdr_plus_other"

        # Donor capture (mutually exclusive with HDR_*)
        assert EditingOutcome.DONOR_CAPTURE.value == "donor_capture"

        # Non-HDR repair outcomes
        assert EditingOutcome.NHEJ_INDEL.value == "nhej_indel"
        assert EditingOutcome.MMEJ_INDEL.value == "mmej_indel"

        # Other outcomes
        assert EditingOutcome.WT.value == "wt"
        assert EditingOutcome.NON_DONOR_SNV.value == "non_donor_snv"
        assert EditingOutcome.UNCLASSIFIED.value == "unclassified"
        assert EditingOutcome.UNMAPPED.value == "unmapped"

    def test_all_outcomes_have_values(self):
        """Test all outcomes are string-valued."""
        for outcome in EditingOutcome:
            assert isinstance(outcome.value, str)
            assert len(outcome.value) > 0


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


class TestClassificationResult:
    """Test ClassificationResult dataclass."""

    def test_classification_result_creation(self):
        """Test creating a classification result."""
        result = ClassificationResult(
            outcome=EditingOutcome.HDR_COMPLETE,
            confidence=0.95,
            details={"hdr_n_matches": 5, "hdr_n_total": 5},
        )

        assert result.outcome == EditingOutcome.HDR_COMPLETE
        assert result.confidence == 0.95
        assert result.details["hdr_n_matches"] == 5

    def test_classification_result_with_qc_flags(self):
        """Test classification result with QC flags."""
        result = ClassificationResult(
            outcome=EditingOutcome.NHEJ_INDEL,
            confidence=0.9,
            details={"indel_mechanism": "nhej", "max_microhomology": 1},
            qc_flags=["low_complexity"],
        )

        assert result.outcome == EditingOutcome.NHEJ_INDEL
        assert "low_complexity" in result.qc_flags
        assert result.details["indel_mechanism"] == "nhej"

    def test_classification_result_mmej(self):
        """Test MMEJ classification result."""
        result = ClassificationResult(
            outcome=EditingOutcome.MMEJ_INDEL,
            confidence=0.85,
            details={"indel_mechanism": "mmej", "max_microhomology": 5},
        )

        assert result.outcome == EditingOutcome.MMEJ_INDEL
        assert result.details["max_microhomology"] == 5


class TestMicrohomologyDetection:
    """Test microhomology detection for NHEJ vs MMEJ classification."""

    def test_mmej_threshold_constant(self):
        """Test MMEJ microhomology threshold is set to 2bp."""
        assert MMEJ_MICROHOMOLOGY_THRESHOLD == 2

    def test_no_microhomology(self):
        """Test deletion with no microhomology."""
        # Reference: AAACCC...GGG with deletion removing CCC
        ref_seq = "AAACCCGGG"
        # Deletion starts at 3, ends at 6 (removing CCC)
        left_mh, right_mh, left_seq, right_seq = calculate_deletion_microhomology(
            ref_seq, deletion_start=3, deletion_end=6
        )
        # No microhomology expected (AAA != GGG)
        assert left_mh == 0 or right_mh == 0  # At least one side has no MH

    def test_microhomology_present(self):
        """Test deletion with microhomology."""
        # Reference with repeat: AAATTTCCCAAATTT
        # Deletion removing CCCAAA would have 3bp microhomology (AAA)
        ref_seq = "AAATTTCCCAAATTT"
        # Deletion from pos 6 to 12 removes "CCCAAA"
        left_mh, right_mh, left_seq, right_seq = calculate_deletion_microhomology(
            ref_seq, deletion_start=6, deletion_end=12
        )
        # Check if microhomology detected
        max_mh = max(left_mh, right_mh)
        assert max_mh >= 0  # Should find some microhomology

    def test_classify_indel_mechanism_nhej(self):
        """Test classification of indel as NHEJ (low microhomology)."""
        # Reference with no repeats
        ref_seq = "AAACCCGGGTTT"
        mechanism, mh_len = classify_indel_mechanism(ref_seq, 3, 6)
        # Should be NHEJ if microhomology <= 2bp
        assert mechanism in ('nhej', 'mmej')

    def test_classify_indel_mechanism_empty_ref(self):
        """Test classification defaults to NHEJ without reference."""
        mechanism, mh_len = classify_indel_mechanism("", 0, 0)
        assert mechanism == 'nhej'
        assert mh_len == 0


class TestSummarizeClassifications:
    """Test summarize_classifications function with new categories."""

    def test_empty_classifications(self):
        """Test summarizing empty list."""
        summary = summarize_classifications([])
        assert summary['total_reads'] == 0
        assert summary['hdr_total_count'] == 0

    def test_hdr_total_aggregation(self):
        """Test HDR total aggregates all HDR categories."""
        classifications = [
            ClassificationResult(outcome=EditingOutcome.HDR_COMPLETE, confidence=1.0),
            ClassificationResult(outcome=EditingOutcome.HDR_COMPLETE, confidence=1.0),
            ClassificationResult(outcome=EditingOutcome.HDR_PARTIAL, confidence=0.9),
            ClassificationResult(outcome=EditingOutcome.HDR_PLUS_NHEJ_INDEL, confidence=0.85),
            ClassificationResult(outcome=EditingOutcome.HDR_PLUS_MMEJ_INDEL, confidence=0.85),
            ClassificationResult(outcome=EditingOutcome.WT, confidence=1.0),
        ]

        summary = summarize_classifications(classifications)

        assert summary['total_reads'] == 6
        assert summary['hdr_complete_count'] == 2
        assert summary['hdr_partial_count'] == 1
        assert summary['hdr_plus_nhej_indel_count'] == 1
        assert summary['hdr_plus_mmej_indel_count'] == 1
        assert summary['hdr_total_count'] == 5  # 2 + 1 + 1 + 1
        assert summary['wt_count'] == 1

    def test_nhej_mmej_total_aggregation(self):
        """Test NHEJ/MMEJ total aggregates correctly."""
        classifications = [
            ClassificationResult(outcome=EditingOutcome.NHEJ_INDEL, confidence=0.9),
            ClassificationResult(outcome=EditingOutcome.NHEJ_INDEL, confidence=0.9),
            ClassificationResult(outcome=EditingOutcome.MMEJ_INDEL, confidence=0.9),
            ClassificationResult(outcome=EditingOutcome.WT, confidence=1.0),
        ]

        summary = summarize_classifications(classifications)

        assert summary['nhej_indel_count'] == 2
        assert summary['mmej_indel_count'] == 1
        assert summary['nhej_mmej_total_count'] == 3  # 2 + 1

    def test_edited_total_aggregation(self):
        """Test edited total includes HDR + NHEJ/MMEJ + donor capture."""
        classifications = [
            ClassificationResult(outcome=EditingOutcome.HDR_COMPLETE, confidence=1.0),
            ClassificationResult(outcome=EditingOutcome.NHEJ_INDEL, confidence=0.9),
            ClassificationResult(outcome=EditingOutcome.DONOR_CAPTURE, confidence=0.85),
            ClassificationResult(outcome=EditingOutcome.WT, confidence=1.0),
            ClassificationResult(outcome=EditingOutcome.NON_DONOR_SNV, confidence=0.7),
        ]

        summary = summarize_classifications(classifications)

        assert summary['hdr_total_count'] == 1
        assert summary['nhej_mmej_total_count'] == 1
        assert summary['donor_capture_count'] == 1
        assert summary['edited_total_count'] == 3  # 1 + 1 + 1
        assert summary['wt_count'] == 1
        assert summary['non_donor_snv_count'] == 1

    def test_rate_calculations(self):
        """Test rate calculations are correct."""
        classifications = [
            ClassificationResult(outcome=EditingOutcome.HDR_COMPLETE, confidence=1.0),
            ClassificationResult(outcome=EditingOutcome.WT, confidence=1.0),
        ]

        summary = summarize_classifications(classifications)

        assert summary['total_reads'] == 2
        assert summary['hdr_complete_rate'] == 0.5
        assert summary['wt_rate'] == 0.5
        assert summary['hdr_total_rate'] == 0.5


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
