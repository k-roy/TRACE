"""
Configuration classes and automatic inference logic for CRISPRo.

Author: Kevin R. Roy
"""

from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import List, Dict, Any, Optional, Tuple
import yaml

from .utils.sequence import reverse_complement, find_guide_in_sequence


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
    """Information about a single edit position."""
    position: int
    ref_base: str
    hdr_base: str

    def __str__(self) -> str:
        return f"Position {self.position}: {self.ref_base} → {self.hdr_base}"


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

    def analyze(self) -> 'LocusConfig':
        """Analyze the locus and populate inferred values."""
        self._find_guide_and_pam()
        self._detect_edits()
        self._detect_homology_arms()
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
        """Detect all edit positions by comparing reference and HDR template."""
        self.edits = []

        # Align HDR template within reference if different lengths
        if len(self.hdr_template) < len(self.reference):
            # Find where HDR template aligns in reference
            hdr_upper = self.hdr_template.upper()
            ref_upper = self.reference.upper()

            # Try to find a matching prefix
            best_offset = 0
            best_matches = 0
            for offset in range(len(self.reference) - len(self.hdr_template) + 1):
                matches = sum(
                    1 for i, base in enumerate(hdr_upper)
                    if ref_upper[offset + i] == base
                )
                if matches > best_matches:
                    best_matches = matches
                    best_offset = offset

            # Compare at best offset
            for i, (ref_base, hdr_base) in enumerate(
                zip(ref_upper[best_offset:best_offset + len(hdr_upper)], hdr_upper)
            ):
                if ref_base != hdr_base:
                    self.edits.append(EditInfo(
                        position=best_offset + i,
                        ref_base=ref_base,
                        hdr_base=hdr_base,
                    ))
        else:
            # Same length, direct comparison
            ref_upper = self.reference.upper()
            hdr_upper = self.hdr_template.upper()
            for i, (ref_base, hdr_base) in enumerate(zip(ref_upper, hdr_upper)):
                if ref_base != hdr_base:
                    self.edits.append(EditInfo(
                        position=i,
                        ref_base=ref_base,
                        hdr_base=hdr_base,
                    ))

    def _detect_homology_arms(self):
        """Detect homology arm boundaries from edit positions."""
        if not self.edits:
            # No edits, entire sequence is homology
            self.homology_arms = HomologyArms(
                left_start=0,
                left_end=len(self.reference),
                left_length=len(self.reference),
                right_start=0,
                right_end=len(self.reference),
                right_length=len(self.reference),
            )
            return

        first_edit = min(e.position for e in self.edits)
        last_edit = max(e.position for e in self.edits)

        self.homology_arms = HomologyArms(
            left_start=0,
            left_end=first_edit,
            left_length=first_edit,
            right_start=last_edit + 1,
            right_end=len(self.reference),
            right_length=len(self.reference) - last_edit - 1,
        )

    def print_summary(self):
        """Print a human-readable summary of the locus configuration."""
        print("\n" + "=" * 60)
        print("=== CRISPRo Analysis Configuration ===")
        print("=" * 60)

        print(f"\nReference sequence: {len(self.reference)} bp")
        print(f"HDR template: {len(self.hdr_template)} bp")

        if self.homology_arms:
            print("\nDonor template analysis:")
            print(f"  - Left homology arm: positions {self.homology_arms.left_start + 1}-"
                  f"{self.homology_arms.left_end} on reference "
                  f"({self.homology_arms.left_length} bp)")
            print(f"  - Right homology arm: positions {self.homology_arms.right_start + 1}-"
                  f"{self.homology_arms.right_end} on reference "
                  f"({self.homology_arms.right_length} bp)")

        if self.edits:
            positions = [str(e.position + 1) for e in self.edits]  # 1-indexed for display
            print(f"  - Donor edits detected at positions: {', '.join(positions)} on reference")
            for edit in self.edits:
                print(f"    * Position {edit.position + 1}: {edit.ref_base} → {edit.hdr_base}")

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
