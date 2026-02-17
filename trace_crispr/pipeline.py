"""
Main pipeline orchestration for CRISPRo.

Author: Kevin R. Roy
"""

import tempfile
from pathlib import Path
from typing import List, Dict, Optional, Set
from dataclasses import dataclass
from concurrent.futures import ProcessPoolExecutor, as_completed
import logging
import pysam

from .config import LocusConfig, PipelineConfig, MultiTemplateLocusConfig
from .io.sample_key import Sample
from .io.output import (
    SampleResult, write_results_tsv, generate_summary_report,
    MultiTemplateSampleResult, write_multi_template_results
)
from .preprocessing.detection import run_auto_detection, AutoDetectionResults
from .preprocessing.trimming import trim_adapters
from .preprocessing.contamination import filter_contamination_fastq, create_contamination_filter
from .integrations.aligners import (
    run_triple_alignment, create_reference_fasta, AlignerManager
)
from .integrations.crispresso import CRISPRessoRunner
from .core.scoring import (
    score_alignment, select_best_alignment, select_best_alignment_paired,
    get_dedup_signature, DeduplicationSignature
)
from .core.classification import (
    classify_read, get_hdr_signature_positions, summarize_classifications,
    EditingOutcome, ClassificationResult
)
from .core.kmer import (
    KmerClassifier, classify_fastq_kmer,
    MultiTemplateKmerClassifier, classify_fastq_multi_template
)

logger = logging.getLogger(__name__)


@dataclass
class PipelineState:
    """Track pipeline state for a sample."""
    sample: Sample
    trimmed_r1: Optional[Path] = None
    trimmed_r2: Optional[Path] = None
    cleaned_r1: Optional[Path] = None
    cleaned_r2: Optional[Path] = None
    alignment_dir: Optional[Path] = None
    best_bam: Optional[Path] = None


class EditingPipeline:
    """Main pipeline orchestrator."""

    def __init__(self, config: PipelineConfig):
        self.config = config
        self.locus = config.locus
        self.aligner_manager = AlignerManager()

        # Prepare HDR signature for classification
        self._prepare_hdr_signature()

        # Prepare k-mer classifier
        self._prepare_kmer_classifier()

    def _prepare_hdr_signature(self):
        """Prepare HDR signature positions for classification."""
        if not self.locus.edits:
            self.hdr_signature = []
            return

        # Create aligned reference and HDR sequences
        ref = self.locus.reference
        hdr = self.locus.hdr_template

        # Ensure same length for signature detection
        min_len = min(len(ref), len(hdr))
        self.hdr_signature = get_hdr_signature_positions(
            ref[:min_len], hdr[:min_len]
        )

        logger.info(f"HDR signature: {len(self.hdr_signature)} positions")

    def _prepare_kmer_classifier(self):
        """Prepare k-mer classifier."""
        edit_positions = [e.position for e in self.locus.edits]

        contamination_seq = self.config.contaminant_sequence

        self.kmer_classifier = KmerClassifier.from_sequences(
            reference=self.locus.reference,
            hdr_template=self.locus.hdr_template,
            edit_positions=edit_positions,
            contamination_sequence=contamination_seq,
            kmer_size=self.config.contaminant_kmer_size,
        )

    def run(self, samples: List[Sample] = None) -> List[SampleResult]:
        """
        Run the full pipeline.

        Args:
            samples: List of samples to process (uses config.samples if None)

        Returns:
            List of SampleResult objects
        """
        if samples is None:
            samples = self.config.samples

        # Check aligner availability
        missing = self.aligner_manager.check_requirements()
        if missing:
            logger.warning(
                f"Missing aligners: {missing}. "
                "Install with: conda install -c bioconda bwa bbmap minimap2"
            )

        # Create output directory
        self.config.output_dir.mkdir(parents=True, exist_ok=True)

        # Create reference FASTA
        ref_fasta = self.config.output_dir / "reference.fasta"
        create_reference_fasta(self.locus.reference, ref_fasta)

        # Process samples
        results = []
        for i, sample in enumerate(samples):
            logger.info(f"Processing sample {i+1}/{len(samples)}: {sample.sample_id}")

            try:
                result = self._process_sample(sample, ref_fasta)
                results.append(result)
            except Exception as e:
                logger.error(f"Error processing {sample.sample_id}: {e}")
                results.append(SampleResult(
                    sample_id=sample.sample_id,
                    metadata=sample.metadata,
                ))

        # Write output
        output_file = self.config.output_dir / "per_sample_editing_outcomes_all_methods.tsv"
        write_results_tsv(results, output_file)

        # Generate summary report
        summary_file = self.config.output_dir / "summary_report.md"
        generate_summary_report(results, summary_file)

        return results

    def _process_sample(
        self,
        sample: Sample,
        ref_fasta: Path,
    ) -> SampleResult:
        """Process a single sample through the pipeline."""
        sample_dir = self.config.output_dir / sample.sample_id
        sample_dir.mkdir(parents=True, exist_ok=True)

        state = PipelineState(sample=sample)

        # Step 1: K-mer classification (pre-alignment)
        logger.info(f"  Running k-mer classification...")
        kmer_results = classify_fastq_kmer(
            sample.r1_path,
            sample.r2_path,
            self.kmer_classifier,
        )

        # Step 2: Adapter trimming
        logger.info(f"  Trimming adapters...")
        trim_dir = sample_dir / "trimmed"
        trim_result = trim_adapters(
            sample.r1_path,
            sample.r2_path,
            trim_dir,
            library_type='auto',
            threads=self.config.threads,
        )

        if trim_result.success:
            state.trimmed_r1 = trim_result.output_r1
            state.trimmed_r2 = trim_result.output_r2
        else:
            logger.warning(f"  Trimming failed: {trim_result.error_message}")
            state.trimmed_r1 = sample.r1_path
            state.trimmed_r2 = sample.r2_path

        # Step 3: Contamination filtering (if configured)
        if self.config.contaminant_sequence:
            logger.info(f"  Filtering contamination...")
            clean_dir = sample_dir / "cleaned"
            clean_dir.mkdir(exist_ok=True)

            contamination_kmers = create_contamination_filter(
                Path(self.config.contaminant_sequence) if isinstance(
                    self.config.contaminant_sequence, str
                ) else None,
                kmer_size=self.config.contaminant_kmer_size,
            ) if isinstance(self.config.contaminant_sequence, str) else set()

            # For now, skip if not a path
            state.cleaned_r1 = state.trimmed_r1
            state.cleaned_r2 = state.trimmed_r2
        else:
            state.cleaned_r1 = state.trimmed_r1
            state.cleaned_r2 = state.trimmed_r2

        # Step 4: Triple alignment
        logger.info(f"  Running triple alignment...")
        align_dir = sample_dir / "alignments"
        alignment_results = run_triple_alignment(
            state.cleaned_r1,
            state.cleaned_r2,
            ref_fasta,
            align_dir,
            threads=self.config.threads,
            aligners=self.config.aligners,
        )

        state.alignment_dir = align_dir

        # Step 5: Select best alignments and classify
        logger.info(f"  Classifying reads...")
        classification_result = self._classify_from_alignments(
            alignment_results,
            sample_dir,
        )

        # Step 6: CRISPResso (if enabled)
        crispresso_hdr = None
        crispresso_nhej = None
        if self.config.run_crispresso:
            logger.info(f"  Running CRISPResso2...")
            crispresso_result = self._run_crispresso(sample, sample_dir)
            if crispresso_result:
                crispresso_hdr = crispresso_result.hdr_rate
                crispresso_nhej = crispresso_result.nhej_rate

        # Build result
        return SampleResult(
            sample_id=sample.sample_id,
            total_reads=kmer_results.total_reads,
            aligned_reads=classification_result.get('aligned_reads', 0),
            classifiable_reads=classification_result.get('classifiable_reads', 0),
            duplicate_rate=classification_result.get('duplicate_rate', 0),
            dedup_wt_count=classification_result.get('wt_count', 0),
            dedup_hdr_count=classification_result.get('hdr_count', 0),
            dedup_nhej_count=classification_result.get('nhej_count', 0),
            dedup_large_del_count=classification_result.get('large_del_count', 0),
            nondedup_hdr_count=classification_result.get('nondedup_hdr_count', 0),
            nondedup_nhej_count=classification_result.get('nondedup_nhej_count', 0),
            kmer_wt_count=kmer_results.wt_count,
            kmer_hdr_count=kmer_results.hdr_count,
            kmer_contamination_count=kmer_results.contamination_count,
            crispresso_hdr_rate=crispresso_hdr,
            crispresso_nhej_rate=crispresso_nhej,
            metadata=sample.metadata,
        )

    def _classify_from_alignments(
        self,
        alignment_results: Dict,
        sample_dir: Path,
    ) -> Dict:
        """Classify reads from alignment results."""
        # Find successful alignments
        successful = {k: v for k, v in alignment_results.items() if v.success}

        if not successful:
            logger.warning("No successful alignments")
            return {}

        # Get cut site from locus config
        cut_site = self.locus.guide_info.cleavage_site if self.locus.guide_info else 0

        # Process reads from all aligners and select best
        classifications = []
        seen_signatures: Set[DeduplicationSignature] = set()
        aligned_reads = 0
        duplicates = 0
        nondedup_hdr = 0
        nondedup_nhej = 0

        # For simplicity, use the first successful aligner
        # Full implementation would iterate all and select best per-read
        aligner_name = next(iter(successful))
        bam_path = successful[aligner_name].bam_path

        try:
            with pysam.AlignmentFile(str(bam_path), "rb") as bam:
                for read in bam.fetch(until_eof=True):
                    if read.is_unmapped:
                        continue

                    aligned_reads += 1

                    # Score and classify
                    score = score_alignment(read, aligner_name)

                    # Classify read
                    result = classify_read(
                        read,
                        self.hdr_signature,
                        cut_site,
                        window_size=self.config.analysis_window,
                        large_del_min_size=self.config.large_deletion_min,
                        hdr_threshold=self.config.hdr_threshold,
                    )

                    # Track non-dedup counts
                    if result.outcome in (EditingOutcome.HDR_PERFECT, EditingOutcome.HDR_IMPERFECT):
                        nondedup_hdr += 1
                    elif result.outcome in (EditingOutcome.NHEJ_INSERTION, EditingOutcome.NHEJ_DELETION):
                        nondedup_nhej += 1

                    # Deduplication (simplified for single-end)
                    dedup_sig = DeduplicationSignature(
                        r1_ref_start=score.ref_start,
                        r1_ref_end=score.ref_end,
                        r2_ref_start=0,
                        r2_ref_end=0,
                        r1_cigar=score.cigar,
                        r2_cigar="",
                    )

                    if dedup_sig in seen_signatures:
                        duplicates += 1
                        continue

                    seen_signatures.add(dedup_sig)
                    classifications.append(result)

        except Exception as e:
            logger.error(f"Error reading BAM: {e}")
            return {}

        # Summarize
        summary = summarize_classifications(classifications)

        unique_reads = len(classifications)
        duplicate_rate = duplicates / (aligned_reads) if aligned_reads > 0 else 0

        return {
            'aligned_reads': aligned_reads,
            'classifiable_reads': unique_reads,
            'duplicate_rate': duplicate_rate,
            'wt_count': summary.get('wild_type_count', 0),
            'hdr_count': summary.get('hdr_total_count', 0),
            'nhej_count': summary.get('nhej_insertion_count', 0) + summary.get('nhej_deletion_count', 0),
            'large_del_count': summary.get('large_deletion_count', 0),
            'nondedup_hdr_count': nondedup_hdr,
            'nondedup_nhej_count': nondedup_nhej,
        }

    def _run_crispresso(self, sample: Sample, sample_dir: Path):
        """Run CRISPResso2 for a sample."""
        runner = CRISPRessoRunner()

        if not runner.is_available():
            logger.warning("CRISPResso2 not available, skipping")
            return None

        crispresso_dir = sample_dir / "crispresso"

        return runner.run(
            r1_fastq=sample.r1_path,
            r2_fastq=sample.r2_path,
            amplicon_seq=self.locus.reference,
            guide_seq=self.locus.guide,
            output_dir=crispresso_dir,
            hdr_seq=self.locus.hdr_template,
            name=sample.sample_id,
        )


def run_pipeline(
    reference: str,
    hdr_template: str,
    guide: str,
    samples: List[Sample],
    output_dir: Path,
    nuclease: str = 'cas9',
    threads: int = 4,
    run_crispresso: bool = True,
    contaminant: Optional[str] = None,
) -> List[SampleResult]:
    """
    Convenience function to run the full pipeline.

    Args:
        reference: Reference sequence
        hdr_template: HDR template sequence
        guide: Guide sequence
        samples: List of Sample objects
        output_dir: Output directory
        nuclease: 'cas9' or 'cas12a'
        threads: Number of threads
        run_crispresso: Run CRISPResso2 validation
        contaminant: Optional contamination sequence

    Returns:
        List of SampleResult objects
    """
    from .config import LocusConfig, NucleaseType

    nuclease_type = NucleaseType.CAS9 if nuclease == 'cas9' else NucleaseType.CAS12A

    locus = LocusConfig(
        name="analysis",
        reference=reference,
        hdr_template=hdr_template,
        guide=guide,
        nuclease=nuclease_type,
    ).analyze()

    config = PipelineConfig(
        locus=locus,
        samples=samples,
        output_dir=output_dir,
        threads=threads,
        run_crispresso=run_crispresso,
        contaminant_sequence=contaminant,
    )

    pipeline = EditingPipeline(config)
    return pipeline.run()


@dataclass
class MultiTemplatePipelineConfig:
    """Configuration for multi-template pipeline."""
    locus: MultiTemplateLocusConfig
    samples: List[Sample]
    output_dir: Path

    # Column in sample metadata containing expected template ID
    expected_template_column: str = "expected_barcode"

    # Library type (auto-detected if None)
    library_type: Optional[str] = None

    # Contaminant filtering
    contaminant_sequence: Optional[str] = None
    contaminant_kmer_size: int = 12

    # Processing options
    threads: int = 4
    kmer_size: int = 12
    min_read_length: int = 50

    # Alignment options
    aligners: List[str] = None

    # Classification thresholds
    hdr_threshold: float = 0.8
    large_deletion_min: int = 21
    analysis_window: int = 20

    # CRISPResso integration
    run_crispresso: bool = True

    def __post_init__(self):
        if self.aligners is None:
            self.aligners = ['bwa', 'bbmap', 'minimap2']


class MultiTemplateEditingPipeline:
    """
    Pipeline supporting multiple HDR templates.

    Used for barcode screening experiments where each sample may have
    a different expected barcode, but we want to check for all possible
    barcodes in each sample.
    """

    def __init__(self, config: MultiTemplatePipelineConfig):
        self.config = config
        self.locus = config.locus
        self.aligner_manager = AlignerManager()

        # Prepare multi-template k-mer classifier
        self._prepare_multi_template_classifier()

    def _prepare_multi_template_classifier(self):
        """Prepare k-mer classifier for all templates."""
        edit_positions_by_template = self.locus.get_edit_positions_by_template()

        # Determine optimal k-mer size
        kmer_size = max(self.config.kmer_size, self.locus.recommended_kmer_size())
        if kmer_size != self.config.kmer_size:
            logger.info(f"Adjusting k-mer size to {kmer_size} based on edit sizes")

        self.kmer_classifier = MultiTemplateKmerClassifier.from_multi_template(
            reference=self.locus.reference,
            hdr_templates=self.locus.hdr_templates,
            edit_positions_by_template=edit_positions_by_template,
            contamination_sequence=self.config.contaminant_sequence,
            kmer_size=kmer_size,
        )

        logger.info(f"Created multi-template classifier with {len(self.locus.hdr_templates)} templates")

    def run(self, samples: List[Sample] = None) -> List[MultiTemplateSampleResult]:
        """
        Run the full multi-template pipeline.

        Args:
            samples: List of samples to process (uses config.samples if None)

        Returns:
            List of MultiTemplateSampleResult objects
        """
        if samples is None:
            samples = self.config.samples

        # Check aligner availability
        missing = self.aligner_manager.check_requirements()
        if missing:
            logger.warning(
                f"Missing aligners: {missing}. "
                "Install with: conda install -c bioconda bwa bbmap minimap2"
            )

        # Create output directory
        self.config.output_dir.mkdir(parents=True, exist_ok=True)

        # Create reference FASTA
        ref_fasta = self.config.output_dir / "reference.fasta"
        create_reference_fasta(self.locus.reference, ref_fasta)

        # Process samples
        results = []
        for i, sample in enumerate(samples):
            logger.info(f"Processing sample {i+1}/{len(samples)}: {sample.sample_id}")

            try:
                result = self._process_sample(sample, ref_fasta)
                results.append(result)
            except Exception as e:
                logger.error(f"Error processing {sample.sample_id}: {e}")
                import traceback
                traceback.print_exc()
                results.append(MultiTemplateSampleResult(
                    sample_id=sample.sample_id,
                    total_reads=0,
                    aligned_reads=0,
                    metadata=sample.metadata,
                ))

        # Write output files
        write_multi_template_results(results, self.config.output_dir)

        logger.info(f"Pipeline complete. Results in {self.config.output_dir}")

        return results

    def _process_sample(
        self,
        sample: Sample,
        ref_fasta: Path,
    ) -> MultiTemplateSampleResult:
        """Process a single sample through the multi-template pipeline."""
        sample_dir = self.config.output_dir / sample.sample_id
        sample_dir.mkdir(parents=True, exist_ok=True)

        # Step 1: Multi-template k-mer classification (pre-alignment)
        logger.info(f"  Running multi-template k-mer classification...")
        kmer_results = classify_fastq_multi_template(
            sample.r1_path,
            sample.r2_path,
            self.kmer_classifier,
        )

        # Get expected template for this sample
        expected_template = sample.metadata.get(
            self.config.expected_template_column
        )
        if expected_template and str(expected_template).upper() in ('NA', 'NAN', ''):
            expected_template = None

        # Build result from k-mer classification
        # (Full implementation would also do alignment-based classification)
        result = MultiTemplateSampleResult(
            sample_id=sample.sample_id,
            total_reads=kmer_results.total_reads,
            aligned_reads=kmer_results.total_reads,  # Approximation for k-mer-only
            hdr_counts_by_template=kmer_results.hdr_counts_by_template.copy(),
            wt_count=kmer_results.wt_count,
            nhej_count=0,  # K-mer classification doesn't detect NHEJ
            large_del_count=0,  # K-mer classification doesn't detect large deletions
            ambiguous_count=kmer_results.ambiguous_count,
            expected_template=str(expected_template) if expected_template else None,
            metadata=sample.metadata,
        )

        return result


def run_multi_template_pipeline(
    reference: str,
    hdr_templates: Dict[str, str],
    guide: str,
    samples: List[Sample],
    output_dir: Path,
    nuclease: str = 'cas9',
    expected_template_column: str = 'expected_barcode',
    threads: int = 4,
    run_crispresso: bool = False,
    contaminant: Optional[str] = None,
) -> List[MultiTemplateSampleResult]:
    """
    Convenience function to run the multi-template pipeline.

    Args:
        reference: Reference sequence
        hdr_templates: Dict mapping template_id → HDR template sequence
        guide: Guide sequence
        samples: List of Sample objects
        output_dir: Output directory
        nuclease: 'cas9' or 'cas12a'
        expected_template_column: Column in sample metadata with expected template
        threads: Number of threads
        run_crispresso: Run CRISPResso2 validation (not yet implemented for multi-template)
        contaminant: Optional contamination sequence

    Returns:
        List of MultiTemplateSampleResult objects
    """
    from .config import MultiTemplateLocusConfig, NucleaseType

    nuclease_type = NucleaseType.CAS9 if nuclease == 'cas9' else NucleaseType.CAS12A

    locus = MultiTemplateLocusConfig(
        name="analysis",
        reference=reference,
        hdr_templates=hdr_templates,
        guide=guide,
        nuclease=nuclease_type,
    ).analyze()

    config = MultiTemplatePipelineConfig(
        locus=locus,
        samples=samples,
        output_dir=output_dir,
        expected_template_column=expected_template_column,
        threads=threads,
        run_crispresso=run_crispresso,
        contaminant_sequence=contaminant,
    )

    pipeline = MultiTemplateEditingPipeline(config)
    return pipeline.run()
