# Changelog

All notable changes to TRACE will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.6.1] - 2026-03-18

### Fixed

- Package metadata and documentation updates from audit.

## [0.6.0] - 2026-03-18

### Fixed

- **Homology arm calculation**: Fixed off-by-one error in right homology arm length calculation. Now correctly accounts for multi-base edits (deletions, insertions) when determining the end of the edit region.

- **UMI deduplication quality bias**: Changed from sum-based quality scoring to mean quality per base. Previously, longer reads were favored over shorter high-quality reads during UMI deduplication.

- **R1/R2 read count validation**: Added explicit validation that paired-end FASTQ files have matching read counts. Mismatched files now raise a clear error instead of silently processing corrupted data.

- **Batch processing silent failures**: Preprocessing failures are now explicitly logged with warnings, and affected samples are flagged with `preprocessing_failed=True` in the output. A summary of failed samples is logged at the end of batch processing.

- **Malformed FASTQ handling**: Added validation for FASTQ format (4 lines per record, proper `+` separator). Malformed files now raise descriptive errors.

- **Empty sequence validation**: Needleman-Wunsch and Smith-Waterman alignment functions now raise `ValueError` for empty input sequences instead of returning empty alignments.

### Added

- **FDR correction `use_total_tests` parameter**: Optional conservative FDR correction that uses total tests (including NaN) instead of only valid tests. Default behavior unchanged.

- **`strict` mode for sample key loading**: New `strict=True` parameter in `load_sample_key()` raises exceptions for validation errors (missing files, etc.) instead of just logging warnings.

- **Support for USEARCH/VSEARCH count format**: Pre-collapsed FASTQs with `;size=N` headers (USEARCH format) are now supported in addition to `;count=N` (TRACE format).

- **`preprocessing_failed` and `preprocessing_error` fields**: `CollapsedSample` dataclass now tracks preprocessing failures for downstream QC.

### Changed

- **README overhaul**: Comprehensive documentation update reflecting all 12 classification categories, NHEJ/MMEJ distinction, output columns, and CLI commands. Removed outdated "MIXED" category documentation.

### Documentation

- Added complete classification category reference with biological descriptions
- Documented NHEJ vs MMEJ microhomology-based classification
- Added full output column reference table
- Added CLI commands reference table

## [0.4.0] - 2026-02-20

### Added

- **Multi-reference alignment classifier** (`trace_crispr.core.multi_ref_classifier`): Classify reads by aligning to a FASTA containing WT + all HDR variants. Provides accurate HDR_PERFECT vs HDR_IMPERFECT distinction with configurable mismatch thresholds.

- **Multi-reference FASTA builder** (`trace_crispr.utils.multi_ref_builder`): Build full-length HDR amplicons from short donor templates using homology arm alignment. Automatically detects anchor sequences and expands templates to full amplicon length.

- **Alignment-only classification mode**: Now the default! Pure alignment-based classification using multi-reference FASTA. More accurate than k-mer classification, especially for experiments with many similar barcodes. Set `alignment_only=False` to use legacy k-mer mode.

- **Global sequence deduplication** with `min_global_count` filter: Remove singleton sequences (likely sequencing errors) across all samples before classification. Reduces false positives and speeds up processing.

- **Primary-only alignment flags** for BWA, BBMap, and minimap2: Suppress secondary alignments when aligning to multi-reference FASTA. Ensures each read gets exactly one best alignment per reference.

- **EditingTemplate data model** (`trace_crispr.core.models`): Structured representation of HDR templates with barcode, homology arms, and full sequence.

### Changed

- **BREAKING**: `alignment_only` parameter now defaults to `True`. Use `alignment_only=False` for legacy k-mer mode.
- BBMap subprocess handling: Redirect stderr to /dev/null to prevent pipe buffer deadlock on large verbose output.

### Fixed

- Version mismatch between `__init__.py` (was 0.2.0) and `pyproject.toml` (was 0.3.1).

## [0.3.1] - 2026-02-17

### Added

- 3-pass UMI detection algorithm for improved merge rates on weak-signal libraries (35% → 92.5% improvement).
- Barcode-style template auto-detection and optimized k-mer generation for better discrimination.
- Multi-template batch processing with global sequence deduplication (~100x faster than per-sample processing).
- CRISPResso2 optional validation dependency for result comparison.

## [0.3.0] - 2026-02-10

### Added

- Multi-HDR template support for barcode screening experiments.
- `trace multi-template` command for batch processing.
- `trace generate-templates` and `trace generate-manifest` helper commands.

## [0.2.0] - 2026-01-20

### Added

- Flexible input: DNA sequences directly or FASTA file paths.
- Large edit support (insertions up to 50+ bp).
- K-mer classification with auto-sizing for improved sensitivity.

## [0.1.0] - 2026-01-12

### Added

- Initial release.
- Triple-aligner consensus (BWA-MEM, BBMap, minimap2).
- Automatic PAM/cleavage site detection for Cas9 and Cas12a.
- UMI-based and position-based PCR deduplication.
- Classification of WT, HDR, NHEJ, and large deletion outcomes.
