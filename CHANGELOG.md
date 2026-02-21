# Changelog

All notable changes to TRACE will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
