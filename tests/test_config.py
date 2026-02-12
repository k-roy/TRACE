"""Tests for trace.config module."""

import pytest
from trace_crispr.config import LocusConfig, NucleaseType


# BFP/GFP test sequences
BFP_REFERENCE = (
    "TGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACC"
    "CACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCC"
    "CGAAGGCTACGTCCAGGAGCGCACCAT"
)

GFP_HDR_TEMPLATE = (
    "TGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACG"
    "TACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCC"
    "CGAAGGCTACGTCCAGGAGCGCACCAT"
)

GUIDE_SEQUENCE = "GCTGAAGCACTGCACGCCGT"


class TestLocusConfig:
    """Test LocusConfig class."""

    def test_basic_initialization(self):
        """Test basic LocusConfig creation."""
        locus = LocusConfig(
            name="test",
            reference=BFP_REFERENCE,
            hdr_template=GFP_HDR_TEMPLATE,
            guide=GUIDE_SEQUENCE,
            nuclease=NucleaseType.CAS9,
        )

        assert locus.name == "test"
        assert locus.reference == BFP_REFERENCE
        assert locus.hdr_template == GFP_HDR_TEMPLATE
        assert locus.guide == GUIDE_SEQUENCE
        assert locus.nuclease == NucleaseType.CAS9

    def test_analyze_detects_edits(self):
        """Test that analyze() detects HDR edits."""
        locus = LocusConfig(
            name="test",
            reference=BFP_REFERENCE,
            hdr_template=GFP_HDR_TEMPLATE,
            guide=GUIDE_SEQUENCE,
            nuclease=NucleaseType.CAS9,
        ).analyze()

        assert locus.edits is not None
        assert len(locus.edits) == 2

        # Check edit positions (0-indexed: 70, 71)
        edit_positions = [e.position for e in locus.edits]
        assert 70 in edit_positions
        assert 71 in edit_positions

    def test_analyze_detects_homology_arms(self):
        """Test that analyze() detects homology arms."""
        locus = LocusConfig(
            name="test",
            reference=BFP_REFERENCE,
            hdr_template=GFP_HDR_TEMPLATE,
            guide=GUIDE_SEQUENCE,
            nuclease=NucleaseType.CAS9,
        ).analyze()

        assert locus.homology_arms is not None
        assert locus.homology_arms.left_start == 0
        assert locus.homology_arms.left_end == 70  # First edit at position 70

    def test_analyze_finds_guide_on_minus_strand(self):
        """Test that analyze() finds the guide on the minus strand."""
        locus = LocusConfig(
            name="test",
            reference=BFP_REFERENCE,
            hdr_template=GFP_HDR_TEMPLATE,
            guide=GUIDE_SEQUENCE,
            nuclease=NucleaseType.CAS9,
        ).analyze()

        assert locus.guide_info is not None
        assert locus.guide_info.strand == '-'
        # Guide is on minus strand, so PAM is upstream

    def test_analyze_calculates_cleavage_site(self):
        """Test that analyze() calculates cleavage site for Cas9."""
        locus = LocusConfig(
            name="test",
            reference=BFP_REFERENCE,
            hdr_template=GFP_HDR_TEMPLATE,
            guide=GUIDE_SEQUENCE,
            nuclease=NucleaseType.CAS9,
        ).analyze()

        assert locus.guide_info is not None
        # Cleavage site should be within reasonable distance of guide
        assert locus.guide_info.cleavage_site > 0

    def test_guide_not_found_raises_error(self):
        """Test that analyze() raises error when guide not found."""
        locus = LocusConfig(
            name="test",
            reference=BFP_REFERENCE,
            hdr_template=GFP_HDR_TEMPLATE,
            guide="NNNNNNNNNNNNNNNNNNNN",  # Non-existent guide
            nuclease=NucleaseType.CAS9,
        )

        with pytest.raises(ValueError, match="not found"):
            locus.analyze()


class TestNucleaseType:
    """Test NucleaseType enum."""

    def test_cas9_value(self):
        """Test Cas9 enum value."""
        assert NucleaseType.CAS9.value == "cas9"

    def test_cas12a_value(self):
        """Test Cas12a enum value."""
        assert NucleaseType.CAS12A.value == "cas12a"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
