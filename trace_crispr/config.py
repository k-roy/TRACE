"""
Configuration classes and automatic inference logic for TRACE.

TRACE: Triple-aligner Read Analysis for CRISPR Editing

Author: Kevin R. Roy
"""

from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import List, Dict, Any, Optional, Tuple
import re
import yaml

from .utils.sequence import reverse_complement, find_guide_in_sequence


# Regex to detect if string is pure DNA sequence
DNA_PATTERN = re.compile(r'^[ACGTacgtNn]+$')


def is_dna_sequence(s: str) -> bool:
    """Check if string is a pure DNA sequence (not a file path)."""
    # Must be non-empty and contain only valid DNA characters
    return bool(s) and bool(DNA_PATTERN.match(s)) and len(s) < 1000


def parse_sequence_input(value: str) -> str:
    """
    Parse sequence input - can be either a DNA string or a FASTA file path.

    Args:
        value: Either a DNA sequence string or path to a FASTA file

    Returns:
        The DNA sequence (uppercase)

    Examples:
        >>> parse_sequence_input("ATCGATCG")
        'ATCGATCG'
        >>> parse_sequence_input("reference.fasta")
        'ATCG...'  # contents of file
    """
    value = value.strip()

    # Check if it's a DNA sequence
    if is_dna_sequence(value):
        return value.upper()

    # Otherwise treat as file path
    path = Path(value)
    if not path.exists():
        raise ValueError(f"File not found: {value}")

    return load_fasta(path)


class NucleaseType(Enum):
    """Supported nuclease types with their cleavage geometry."""
    CAS9 = "cas9"
    CAS12A = "cas12a"


@dataclass
class NucleaseConfig:
    """Configuration for a nuclease type."""
    name: str
    pam_pattern: str
    pam_position: str  # '3prime' or '5prime' relative to guide
    cleavage_offset: int  # bp from PAM to primary cleavage site
    cleavage_type: str  # 'blunt' or 'staggered'
    stagger_size: Optional[int] = None  # For staggered cuts

    @classmethod
    def cas9(cls) -> 'NucleaseConfig':
        """SpCas9 configuration."""
        return cls(
            name="SpCas9",
            pam_pattern="NGG",
            pam_position="3prime",
            cleavage_offset=-3,  # 3bp upstream of PAM
            cleavage_type="blunt",
        )

    @classmethod
    def cas12a(cls) -> 'NucleaseConfig':
        """LbCas12a (Cpf1) configuration."""
        return cls(
            name="LbCas12a",
            pam_pattern="TTTN",
            pam_position="5prime",
            cleavage_offset=18,  # 18-19bp downstream of PAM on target strand
            cleavage_type="staggered",
            stagger_size=5,  # 4-5nt overhang
        )


@dataclass
class EditInfo:
    """Information about a single edit at a position."""
    position: int  # Position in reference
    ref_bases: str  # Reference base(s) - can be multiple for deletions or empty for insertions
    hdr_bases: str  # HDR base(s) - can be multiple for insertions or empty for deletions
    edit_type: str = "substitution"  # 'substitution', 'insertion', or 'deletion'

    @property
    def size(self) -> int:
        """Size of the edit (positive for insertions, negative for deletions, 0 for subs)."""
        return len(self.hdr_bases) - len(self.ref_bases)

    def __str__(self) -> str:
        if self.edit_type == "substitution":
            return f"Position {self.position}: {self.ref_bases} → {self.hdr_bases}"
        elif self.edit_type == "insertion":
            return f"Position {self.position}: insert {self.hdr_bases} ({len(self.hdr_bases)} bp)"
        else:  # deletion
            return f"Position {self.position}: delete {self.ref_bases} ({len(self.ref_bases)} bp)"


@dataclass
class GuideInfo:
    """Information about guide targeting."""
    sequence: str
    start: int
    end: int
    strand: str  # '+' or '-'
    pam_start: int
    pam_end: int
    pam_seq: str
    cleavage_site: int


@dataclass
class HomologyArms:
    """Homology arm information from donor template."""
    left_start: int
    left_end: int
    left_length: int
    right_start: int
    right_end: int
    right_length: int


@dataclass
class LocusConfig:
    """Complete locus configuration with inferred values."""
    name: str
    reference: str
    hdr_template: str
    guide: str
    nuclease: NucleaseType

    # Inferred values (populated by analyze())
    guide_info: Optional[GuideInfo] = None
    edits: List[EditInfo] = field(default_factory=list)
    homology_arms: Optional[HomologyArms] = None
    template_offset: int = 0  # Where template aligns in reference

    def analyze(self) -> 'LocusConfig':
        """Analyze the locus and populate inferred values."""
        self._find_guide_and_pam()
        self._detect_edits()
        self._detect_homology_arms()
        return self

    @property
    def max_edit_size(self) -> int:
        """Get the maximum edit size (absolute value of insertion or deletion)."""
        if not self.edits:
            return 0
        return max(abs(e.size) for e in self.edits)

    @property
    def total_edit_span(self) -> int:
        """Get total span of edits including insertions/deletions."""
        if not self.edits:
            return 0
        # For insertions, count the inserted bases; for deletions/subs, count ref bases
        return sum(max(len(e.ref_bases), len(e.hdr_bases)) for e in self.edits)

    def recommended_kmer_size(self, min_kmer: int = 12, buffer: int = 10) -> int:
        """
        Calculate recommended k-mer size based on edit size.

        K-mer size should be at least 10 bp larger than the maximum edit size
        to ensure k-mers can span the edit site unambiguously.

        Args:
            min_kmer: Minimum k-mer size (default: 12)
            buffer: Additional bp beyond max edit size (default: 10)

        Returns:
            Recommended k-mer size
        """
        max_edit = self.max_edit_size
        recommended = max(min_kmer, max_edit + buffer)
        return recommended

    def _find_guide_and_pam(self):
        """Find guide position and infer PAM/cleavage site."""
        guide_pos, strand = find_guide_in_sequence(self.reference, self.guide)

        nuclease_config = (
            NucleaseConfig.cas9() if self.nuclease == NucleaseType.CAS9
            else NucleaseConfig.cas12a()
        )

        guide_len = len(self.guide)

        if nuclease_config.pam_position == "3prime":
            # PAM is downstream of guide (Cas9)
            if strand == '+':
                pam_start = guide_pos + guide_len
                pam_end = pam_start + len(nuclease_config.pam_pattern)
                # Cleavage is 3bp upstream of PAM
                cleavage = pam_start + nuclease_config.cleavage_offset
            else:
                # Guide is on minus strand, PAM is upstream in reference coords
                pam_end = guide_pos
                pam_start = pam_end - len(nuclease_config.pam_pattern)
                cleavage = pam_end - nuclease_config.cleavage_offset
        else:
            # PAM is upstream of guide (Cas12a)
            if strand == '+':
                pam_end = guide_pos
                pam_start = pam_end - len(nuclease_config.pam_pattern)
                # Cleavage is 18bp downstream of PAM
                cleavage = pam_start + nuclease_config.cleavage_offset
            else:
                pam_start = guide_pos + guide_len
                pam_end = pam_start + len(nuclease_config.pam_pattern)
                cleavage = pam_end - nuclease_config.cleavage_offset

        pam_seq = self.reference[pam_start:pam_end].upper()

        self.guide_info = GuideInfo(
            sequence=self.guide,
            start=guide_pos,
            end=guide_pos + guide_len,
            strand=strand,
            pam_start=pam_start,
            pam_end=pam_end,
            pam_seq=pam_seq,
            cleavage_site=cleavage,
        )

    def _detect_edits(self):
        """
        Detect all edit positions by comparing reference and HDR template.

        Handles:
        - Template shorter than reference (typical: 150bp template in 250bp amplicon)
        - Single nucleotide substitutions
        - Insertions (HDR template has extra bases)
        - Deletions (HDR template is missing bases from reference)
        """
        self.edits = []
        ref_upper = self.reference.upper()
        hdr_upper = self.hdr_template.upper()

        # Find where HDR template best aligns in reference
        self.template_offset = self._find_template_alignment(ref_upper, hdr_upper)

        # Calculate where HDR template ends in reference coordinates
        # This is used to properly bound homology arm calculations
        self.template_end_in_ref = self.template_offset + len(hdr_upper)

        # Get the reference region that aligns with the template
        # IMPORTANT: Don't extend beyond template length to avoid spurious deletions
        # The Needleman-Wunsch alignment handles insertions within the template
        ref_region = ref_upper[self.template_offset:self.template_end_in_ref]

        # Align sequences to detect substitutions, insertions, and deletions
        aligned_ref, aligned_hdr = self._align_sequences(ref_region, hdr_upper)

        # Parse alignment to extract edits
        self._parse_alignment_for_edits(aligned_ref, aligned_hdr)

    def _find_template_alignment(self, ref: str, hdr: str) -> int:
        """
        Find the best position to align HDR template within reference.

        Tries multiple anchors (left, middle, right) to find the best alignment.
        This handles cases where the left portion of the template may not be unique.
        """
        max_offset = max(0, len(ref) - len(hdr) + 50)  # Allow some slack for indels

        # Try different anchor positions
        anchor_positions = [
            (0, min(30, len(hdr) // 4)),                           # Left anchor
            (len(hdr) // 4, len(hdr) // 4 + min(30, len(hdr) // 4)),    # Left-center
            (len(hdr) // 2 - 15, len(hdr) // 2 + 15),              # Center anchor
            (max(0, len(hdr) - 30), len(hdr)),                     # Right anchor
        ]

        best_offset = 0
        best_score = 0

        for anchor_start, anchor_end in anchor_positions:
            anchor = hdr[anchor_start:anchor_end]
            if len(anchor) < 10:
                continue

            # Try exact match first
            pos = ref.find(anchor)
            if pos != -1:
                # Found exact match - calculate offset where template would start
                template_start = pos - anchor_start
                if 0 <= template_start <= max_offset:
                    # Verify this alignment by scoring full overlap
                    score = self._score_alignment(ref, hdr, template_start)
                    if score > best_score:
                        best_score = score
                        best_offset = template_start

            # Also try approximate matching for this anchor region
            for offset in range(max_offset + 1):
                ref_pos = offset + anchor_start
                if ref_pos + len(anchor) > len(ref):
                    break
                ref_region = ref[ref_pos:ref_pos + len(anchor)]
                matches = sum(1 for a, b in zip(anchor, ref_region) if a == b)
                if matches > len(anchor) * 0.8:  # At least 80% match
                    score = self._score_alignment(ref, hdr, offset)
                    if score > best_score:
                        best_score = score
                        best_offset = offset

        return best_offset

    def _score_alignment(self, ref: str, hdr: str, offset: int) -> int:
        """Score how well the template aligns at a given offset."""
        compare_len = min(len(hdr), len(ref) - offset)
        matches = sum(
            1 for i in range(compare_len)
            if ref[offset + i] == hdr[i]
        )
        return matches

    def _align_sequences(self, ref: str, hdr: str) -> Tuple[str, str]:
        """
        Simple Needleman-Wunsch alignment for detecting indels.

        Returns aligned sequences with gaps represented as '-'.
        """
        # Simple implementation for small sequences
        # For large indels, we need proper alignment

        # If lengths are similar (within 5bp), do direct comparison
        if abs(len(ref) - len(hdr)) <= 5 and len(ref) < 500:
            return self._simple_align(ref[:len(hdr)], hdr)

        # For larger differences, use dynamic programming
        return self._needleman_wunsch(ref, hdr)

    def _simple_align(self, ref: str, hdr: str) -> Tuple[str, str]:
        """Simple alignment for sequences of similar length."""
        # Pad shorter sequence
        max_len = max(len(ref), len(hdr))
        ref_padded = ref.ljust(max_len, '-')
        hdr_padded = hdr.ljust(max_len, '-')
        return ref_padded, hdr_padded

    def _needleman_wunsch(self, ref: str, hdr: str) -> Tuple[str, str]:
        """
        Needleman-Wunsch global alignment for indel detection.

        Optimized for HDR template alignment where we expect:
        - Long matching homology arms
        - Short edit region in the middle (substitutions, insertions, or deletions)

        Uses affine gap scoring to prefer fewer, longer gaps over many small gaps.
        """
        m, n = len(ref), len(hdr)

        # Limit alignment length to avoid memory issues
        max_len = min(m, len(hdr) + 100)
        ref = ref[:max_len]
        m = len(ref)

        # Scoring - balanced for typical HDR edits
        # Make mismatches costly so alignment prefers gaps for long insertions
        match_score = 3
        mismatch_score = -3  # Costly mismatches
        gap_score = -1       # Cheaper gaps to encourage longer gaps over mismatches

        # Initialize score matrix
        score = [[0] * (n + 1) for _ in range(m + 1)]

        # Initialize first row and column
        for i in range(m + 1):
            score[i][0] = i * gap_score
        for j in range(n + 1):
            score[0][j] = j * gap_score

        # Fill matrix
        for i in range(1, m + 1):
            for j in range(1, n + 1):
                match = score[i-1][j-1] + (match_score if ref[i-1] == hdr[j-1] else mismatch_score)
                delete = score[i-1][j] + gap_score
                insert = score[i][j-1] + gap_score
                score[i][j] = max(match, delete, insert)

        # Traceback
        aligned_ref = []
        aligned_hdr = []
        i, j = m, n

        while i > 0 or j > 0:
            if i > 0 and j > 0:
                current = score[i][j]
                diag = score[i-1][j-1]
                match = match_score if ref[i-1] == hdr[j-1] else mismatch_score

                if current == diag + match:
                    aligned_ref.append(ref[i-1])
                    aligned_hdr.append(hdr[j-1])
                    i -= 1
                    j -= 1
                elif i > 0 and current == score[i-1][j] + gap_score:
                    aligned_ref.append(ref[i-1])
                    aligned_hdr.append('-')
                    i -= 1
                else:
                    aligned_ref.append('-')
                    aligned_hdr.append(hdr[j-1])
                    j -= 1
            elif i > 0:
                aligned_ref.append(ref[i-1])
                aligned_hdr.append('-')
                i -= 1
            else:
                aligned_ref.append('-')
                aligned_hdr.append(hdr[j-1])
                j -= 1

        return ''.join(reversed(aligned_ref)), ''.join(reversed(aligned_hdr))

    def _parse_alignment_for_edits(self, aligned_ref: str, aligned_hdr: str):
        """Parse aligned sequences to extract edit information."""
        ref_pos = self.template_offset  # Position in original reference
        i = 0

        while i < len(aligned_ref):
            ref_base = aligned_ref[i]
            hdr_base = aligned_hdr[i]

            if ref_base == hdr_base:
                # Match - no edit
                if ref_base != '-':
                    ref_pos += 1
                i += 1
            elif ref_base == '-':
                # Insertion in HDR template
                inserted_bases = []
                while i < len(aligned_ref) and aligned_ref[i] == '-':
                    inserted_bases.append(aligned_hdr[i])
                    i += 1
                self.edits.append(EditInfo(
                    position=ref_pos,
                    ref_bases='',
                    hdr_bases=''.join(inserted_bases),
                    edit_type='insertion'
                ))
            elif hdr_base == '-':
                # Deletion in HDR template
                deleted_bases = []
                while i < len(aligned_ref) and aligned_hdr[i] == '-':
                    deleted_bases.append(aligned_ref[i])
                    i += 1
                    ref_pos += 1
                self.edits.append(EditInfo(
                    position=ref_pos - len(deleted_bases),
                    ref_bases=''.join(deleted_bases),
                    hdr_bases='',
                    edit_type='deletion'
                ))
            else:
                # Substitution
                self.edits.append(EditInfo(
                    position=ref_pos,
                    ref_bases=ref_base,
                    hdr_bases=hdr_base,
                    edit_type='substitution'
                ))
                ref_pos += 1
                i += 1

    def _detect_homology_arms(self):
        """Detect homology arm boundaries from edit positions.

        Homology arms are calculated relative to the HDR template's position
        in the reference, not the full reference length. This properly handles
        shorter HDR templates (e.g., 150bp template in 250bp amplicon).
        """
        # Use template bounds, not full reference
        template_start = self.template_offset
        template_end = getattr(self, 'template_end_in_ref',
                               self.template_offset + len(self.hdr_template))

        if not self.edits:
            # No edits, entire template region is homology
            self.homology_arms = HomologyArms(
                left_start=template_start,
                left_end=template_end,
                left_length=template_end - template_start,
                right_start=template_start,
                right_end=template_end,
                right_length=template_end - template_start,
            )
            return

        first_edit = min(e.position for e in self.edits)
        last_edit = max(e.position for e in self.edits)

        self.homology_arms = HomologyArms(
            left_start=template_start,
            left_end=first_edit,
            left_length=first_edit - template_start,
            right_start=last_edit + 1,
            right_end=template_end,  # Use template end, not reference end
            right_length=template_end - last_edit - 1,
        )

    def print_summary(self):
        """Print a human-readable summary of the locus configuration."""
        print("\n" + "=" * 60)
        print("=== TRACE Analysis Configuration ===")
        print("=" * 60)

        print(f"\nReference sequence: {len(self.reference)} bp")
        print(f"HDR template: {len(self.hdr_template)} bp")

        if self.template_offset > 0:
            print(f"  - Template aligns at position {self.template_offset + 1} in reference")

        if self.homology_arms:
            print("\nDonor template analysis:")
            print(f"  - Left homology arm: positions {self.homology_arms.left_start + 1}-"
                  f"{self.homology_arms.left_end} on reference "
                  f"({self.homology_arms.left_length} bp)")
            print(f"  - Right homology arm: positions {self.homology_arms.right_start + 1}-"
                  f"{self.homology_arms.right_end} on reference "
                  f"({self.homology_arms.right_length} bp)")

        if self.edits:
            print(f"\n  Edits detected ({len(self.edits)} total):")
            for edit in self.edits:
                if edit.edit_type == "substitution":
                    print(f"    * Position {edit.position + 1}: {edit.ref_bases} → {edit.hdr_bases} (substitution)")
                elif edit.edit_type == "insertion":
                    print(f"    * Position {edit.position + 1}: +{edit.hdr_bases} ({len(edit.hdr_bases)} bp insertion)")
                else:  # deletion
                    print(f"    * Position {edit.position + 1}: -{edit.ref_bases} ({len(edit.ref_bases)} bp deletion)")

            # Show edit summary
            max_edit = self.max_edit_size
            if max_edit > 0:
                print(f"\n  Maximum edit size: {max_edit} bp")
                print(f"  Recommended k-mer size: {self.recommended_kmer_size()} bp")

        if self.guide_info:
            print("\nGuide analysis:")
            print(f"  - Guide sequence: {self.guide_info.sequence}")
            print(f"  - Guide targets: positions {self.guide_info.start + 1}-"
                  f"{self.guide_info.end} on reference ({self.guide_info.strand} strand)")
            print(f"  - PAM: {self.guide_info.pam_seq} at positions "
                  f"{self.guide_info.pam_start + 1}-{self.guide_info.pam_end} on reference")
            print(f"  - Cleavage site: position {self.guide_info.cleavage_site + 1} on reference")

        print()


@dataclass
class MultiTemplateLocusConfig:
    """
    Locus configuration with multiple possible HDR templates.

    Used for barcode screening experiments where each sample may have
    a different expected barcode, but we want to check for all possible
    barcodes in each sample (purity/sanity check).
    """
    name: str
    reference: str
    hdr_templates: Dict[str, str]  # template_id → sequence
    guide: str
    nuclease: NucleaseType

    # Inferred per-template (populated by analyze())
    template_edits: Dict[str, List[EditInfo]] = field(default_factory=dict)
    template_offsets: Dict[str, int] = field(default_factory=dict)
    template_ends: Dict[str, int] = field(default_factory=dict)
    guide_info: Optional[GuideInfo] = None
    homology_arms: Optional[HomologyArms] = None

    def analyze(self) -> 'MultiTemplateLocusConfig':
        """
        Analyze all templates and populate inferred values.

        For each template, detects:
        - Template offset (where it aligns in reference)
        - Edits (substitutions, insertions, deletions)
        - Template end position

        Also detects guide position and homology arms (shared across templates).
        """
        # First, find guide position (same for all templates)
        self._find_guide_and_pam()

        # Process each template
        for template_id, template_seq in self.hdr_templates.items():
            # Create a temporary LocusConfig to leverage existing analysis
            temp_locus = LocusConfig(
                name=f"{self.name}_{template_id}",
                reference=self.reference,
                hdr_template=template_seq,
                guide=self.guide,
                nuclease=self.nuclease,
            )
            temp_locus.analyze()

            # Extract per-template values
            self.template_edits[template_id] = temp_locus.edits
            self.template_offsets[template_id] = temp_locus.template_offset
            self.template_ends[template_id] = getattr(
                temp_locus, 'template_end_in_ref',
                temp_locus.template_offset + len(template_seq)
            )

            # Use first template's homology arms as representative
            if self.homology_arms is None and temp_locus.homology_arms:
                self.homology_arms = temp_locus.homology_arms

        return self

    def _find_guide_and_pam(self):
        """Find guide position and infer PAM/cleavage site."""
        guide_pos, strand = find_guide_in_sequence(self.reference, self.guide)

        nuclease_config = (
            NucleaseConfig.cas9() if self.nuclease == NucleaseType.CAS9
            else NucleaseConfig.cas12a()
        )

        guide_len = len(self.guide)

        if nuclease_config.pam_position == "3prime":
            if strand == '+':
                pam_start = guide_pos + guide_len
                pam_end = pam_start + len(nuclease_config.pam_pattern)
                cleavage = pam_start + nuclease_config.cleavage_offset
            else:
                pam_end = guide_pos
                pam_start = pam_end - len(nuclease_config.pam_pattern)
                cleavage = pam_end - nuclease_config.cleavage_offset
        else:
            if strand == '+':
                pam_end = guide_pos
                pam_start = pam_end - len(nuclease_config.pam_pattern)
                cleavage = pam_start + nuclease_config.cleavage_offset
            else:
                pam_start = guide_pos + guide_len
                pam_end = pam_start + len(nuclease_config.pam_pattern)
                cleavage = pam_end - nuclease_config.cleavage_offset

        pam_seq = self.reference[pam_start:pam_end].upper()

        self.guide_info = GuideInfo(
            sequence=self.guide,
            start=guide_pos,
            end=guide_pos + guide_len,
            strand=strand,
            pam_start=pam_start,
            pam_end=pam_end,
            pam_seq=pam_seq,
            cleavage_site=cleavage,
        )

    @property
    def max_edit_size(self) -> int:
        """Get the maximum edit size across all templates."""
        max_size = 0
        for edits in self.template_edits.values():
            for edit in edits:
                max_size = max(max_size, abs(edit.size))
        return max_size

    def recommended_kmer_size(self, min_kmer: int = 12, buffer: int = 10) -> int:
        """Calculate recommended k-mer size based on max edit size."""
        return max(min_kmer, self.max_edit_size + buffer)

    def get_edit_positions_by_template(self) -> Dict[str, List[int]]:
        """Get list of edit positions for each template."""
        result = {}
        for template_id, edits in self.template_edits.items():
            result[template_id] = [e.position for e in edits]
        return result

    def print_summary(self):
        """Print a human-readable summary of the multi-template locus."""
        print("\n" + "=" * 60)
        print("=== TRACE Multi-Template Analysis Configuration ===")
        print("=" * 60)

        print(f"\nReference sequence: {len(self.reference)} bp")
        print(f"Number of HDR templates: {len(self.hdr_templates)}")

        # Show template lengths
        template_lengths = [len(t) for t in self.hdr_templates.values()]
        if len(set(template_lengths)) == 1:
            print(f"Template length: {template_lengths[0]} bp (all same)")
        else:
            print(f"Template lengths: {min(template_lengths)}-{max(template_lengths)} bp")

        # Show guide info
        if self.guide_info:
            print("\nGuide analysis:")
            print(f"  - Guide sequence: {self.guide_info.sequence}")
            print(f"  - Guide targets: positions {self.guide_info.start + 1}-"
                  f"{self.guide_info.end} on reference ({self.guide_info.strand} strand)")
            print(f"  - PAM: {self.guide_info.pam_seq} at positions "
                  f"{self.guide_info.pam_start + 1}-{self.guide_info.pam_end}")
            print(f"  - Cleavage site: position {self.guide_info.cleavage_site + 1}")

        # Show edit summary
        if self.template_edits:
            print(f"\nEdits by template (showing first 5):")
            for i, (template_id, edits) in enumerate(self.template_edits.items()):
                if i >= 5:
                    print(f"  ... and {len(self.template_edits) - 5} more templates")
                    break
                edit_summary = ", ".join(
                    f"{e.edit_type}@{e.position + 1}"
                    for e in edits[:3]
                )
                if len(edits) > 3:
                    edit_summary += f" (+{len(edits) - 3} more)"
                print(f"  - {template_id}: {edit_summary}")

            max_edit = self.max_edit_size
            if max_edit > 0:
                print(f"\n  Maximum edit size: {max_edit} bp")
                print(f"  Recommended k-mer size: {self.recommended_kmer_size()} bp")

        print()

    @classmethod
    def from_fasta(
        cls,
        name: str,
        reference: str,
        hdr_templates_fasta: Path,
        guide: str,
        nuclease: NucleaseType,
    ) -> 'MultiTemplateLocusConfig':
        """
        Create MultiTemplateLocusConfig from a FASTA file of templates.

        Args:
            name: Locus name
            reference: Reference sequence
            hdr_templates_fasta: Path to FASTA with HDR templates
            guide: Guide sequence
            nuclease: Nuclease type

        Returns:
            MultiTemplateLocusConfig with templates loaded
        """
        templates = {}
        current_id = None
        current_seq = []

        with open(hdr_templates_fasta) as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_id is not None:
                        templates[current_id] = ''.join(current_seq).upper()
                    current_id = line[1:].split()[0]
                    current_seq = []
                else:
                    current_seq.append(line)

            if current_id is not None:
                templates[current_id] = ''.join(current_seq).upper()

        return cls(
            name=name,
            reference=reference,
            hdr_templates=templates,
            guide=guide,
            nuclease=nuclease,
        )


@dataclass
class SampleConfig:
    """Sample-level configuration."""
    sample_id: str
    r1_path: Path
    r2_path: Optional[Path] = None
    metadata: Dict[str, Any] = field(default_factory=dict)

    @classmethod
    def from_dict(cls, d: Dict) -> 'SampleConfig':
        """Create from dictionary."""
        return cls(
            sample_id=d['sample_id'],
            r1_path=Path(d['r1_path']),
            r2_path=Path(d['r2_path']) if d.get('r2_path') else None,
            metadata={k: v for k, v in d.items()
                     if k not in ('sample_id', 'r1_path', 'r2_path')},
        )


@dataclass
class PipelineConfig:
    """Full pipeline configuration."""
    locus: LocusConfig
    samples: List[SampleConfig]
    output_dir: Path

    # Library type (auto-detected if None)
    library_type: Optional[str] = None  # 'TruSeq' or 'Tn5'

    # Contaminant filtering
    contaminant_sequence: Optional[str] = None
    contaminant_kmer_size: int = 12

    # Processing options
    threads: int = 4
    min_read_length: int = 50

    # Alignment options
    aligners: List[str] = field(default_factory=lambda: ['bwa', 'bbmap', 'minimap2'])

    # Classification thresholds
    hdr_threshold: float = 0.8  # Fraction of HDR positions required
    large_deletion_min: int = 21  # bp
    analysis_window: int = 20  # bp around cut site

    # CRISPResso integration
    run_crispresso: bool = True

    @classmethod
    def from_yaml(cls, path: Path) -> 'PipelineConfig':
        """Load configuration from YAML file."""
        with open(path) as f:
            data = yaml.safe_load(f)

        # Load locus configuration
        if 'locus_file' in data:
            locus_data = yaml.safe_load(open(data['locus_file']))
        else:
            locus_data = data.get('locus', {})

        # Read sequences from files if paths provided
        reference = locus_data.get('reference', '')
        if reference.endswith(('.fa', '.fasta')):
            reference = _read_fasta_sequence(reference)

        hdr_template = locus_data.get('hdr_template', '')
        if hdr_template.endswith(('.fa', '.fasta')):
            hdr_template = _read_fasta_sequence(hdr_template)

        nuclease_str = data.get('nuclease', 'cas9').lower()
        nuclease = NucleaseType.CAS9 if nuclease_str == 'cas9' else NucleaseType.CAS12A

        locus = LocusConfig(
            name=locus_data.get('name', 'unnamed'),
            reference=reference,
            hdr_template=hdr_template,
            guide=locus_data.get('guide', data.get('guide', '')),
            nuclease=nuclease,
        ).analyze()

        # Load samples
        samples = []
        if 'sample_key' in data:
            import pandas as pd
            df = pd.read_csv(data['sample_key'], sep='\t')
            for _, row in df.iterrows():
                samples.append(SampleConfig.from_dict(row.to_dict()))

        # Read contaminant sequence if provided
        contaminant_seq = None
        if data.get('contaminant'):
            contaminant_path = data['contaminant']
            if Path(contaminant_path).exists():
                contaminant_seq = _read_fasta_sequence(contaminant_path)

        return cls(
            locus=locus,
            samples=samples,
            output_dir=Path(data.get('output_dir', './results')),
            library_type=data.get('library_type'),
            contaminant_sequence=contaminant_seq,
            threads=data.get('threads', 4),
            run_crispresso=data.get('crispresso', True),
        )


def _read_fasta_sequence(path: str) -> str:
    """Read first sequence from a FASTA file."""
    sequence = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if sequence:
                    break  # Only read first sequence
                continue
            sequence.append(line.upper())
    return ''.join(sequence)


def load_fasta(path: Path) -> str:
    """Load sequence from FASTA file."""
    return _read_fasta_sequence(str(path))
