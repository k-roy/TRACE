# TRACE

**T**riple-aligner **R**ead **A**nalysis for **C**RISPR **E**diting

Robust quantification of CRISPR editing outcomes from amplicon sequencing data using consensus alignment across BWA-MEM, BBMap, and minimap2.

[![PyPI version](https://badge.fury.io/py/trace-crispr.svg)](https://pypi.org/project/trace-crispr/)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/trace-crispr.svg)](https://bioconda.github.io/recipes/trace-crispr/README.html)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Quick Start

```bash
# Install
pip install trace-crispr

# Run on a single sample
trace run \
  --reference amplicon.fasta \
  --hdr-template donor.fasta \
  --guide GCTGAAGCACTGCACGCCGT \
  --r1 sample_R1.fastq.gz \
  --r2 sample_R2.fastq.gz \
  --output results/

# Run on multiple samples
trace run \
  --reference amplicon.fasta \
  --hdr-template donor.fasta \
  --guide GCTGAAGCACTGCACGCCGT \
  --sample-key samples.tsv \
  --output results/ \
  --threads 16
```

## Output

### Editing Summary

TRACE produces a comprehensive per-sample summary table with editing outcomes:

```
sample     total_reads  HDR_COMPLETE_%  HDR_PARTIAL_%  NHEJ_INDEL_%  MMEJ_INDEL_%  WT_%
──────────────────────────────────────────────────────────────────────────────────────────
sample_01  125,432      42.1            3.1            8.2           4.1           38.1
sample_02  98,765       48.7            3.4            5.6           3.1           35.8
sample_03  112,890      45.2            3.7            9.8           5.4           31.2
```

### Classification Categories

TRACE classifies reads into 12 comprehensive categories:

**HDR Categories (successful donor integration):**
| Outcome | Description |
|---------|-------------|
| `HDR_COMPLETE` | All donor-encoded edits present, no other modifications |
| `HDR_PARTIAL` | Subset of donor SNVs integrated, no other modifications |
| `HDR_PLUS_NHEJ_INDEL` | Donor edits + classical NHEJ indel at cut site (0-2bp microhomology) |
| `HDR_PLUS_MMEJ_INDEL` | Donor edits + MMEJ indel at cut site (>2bp microhomology) |
| `HDR_PLUS_OTHER` | Donor edits + modifications NOT at cut site |

**Non-HDR Repair Outcomes:**
| Outcome | Description |
|---------|-------------|
| `DONOR_CAPTURE` | Donor edits present + extra donor sequence duplicated at site (over-integration) |
| `NHEJ_INDEL` | Classical NHEJ indel at cut site (0-2bp microhomology at deletion boundaries) |
| `MMEJ_INDEL` | Microhomology-mediated indel (>2bp microhomology at deletion boundaries) |

**Other Outcomes:**
| Outcome | Description |
|---------|-------------|
| `WT` | Wild-type / unedited |
| `NON_DONOR_SNV` | SNVs not matching donor template, no indels |
| `UNCLASSIFIED` | Does not fit other categories |
| `UNMAPPED` | Read did not align to reference |

**Note:** `LARGE_DELETION` (≥50bp) is a flag applied to reads within the above categories, not a separate category.

### NHEJ vs MMEJ Classification

TRACE distinguishes between two types of error-prone repair based on deletion microhomology:

- **NHEJ (Non-Homologous End Joining):** Deletions with 0-2bp of microhomology at the deletion boundaries. This is classical error-prone repair.
- **MMEJ (Microhomology-Mediated End Joining):** Deletions with >2bp of microhomology. This alternative repair pathway uses short homologous sequences to rejoin the break.

This distinction is biologically significant as MMEJ and NHEJ involve different repair proteins and have different mutational signatures.

### Output Columns

The main output file (`per_sample_editing_outcomes_all_methods.tsv`) contains:

| Column | Description |
|--------|-------------|
| `sample` | Sample identifier |
| `total_reads` | Total reads processed |
| `aligned_reads` | Reads that aligned to reference |
| `classifiable_reads` | Reads that could be classified |
| `duplicate_rate` | Fraction of reads removed as duplicates |
| `HDR_COMPLETE_%` | Percentage with complete HDR |
| `HDR_PARTIAL_%` | Percentage with partial HDR |
| `HDR_PLUS_NHEJ_%` | HDR + classical NHEJ indel |
| `HDR_PLUS_MMEJ_%` | HDR + microhomology-mediated indel |
| `HDR_PLUS_OTHER_%` | HDR + other modifications |
| `HDR_total_%` | Sum of all HDR categories |
| `DONOR_CAPTURE_%` | Donor over-integration |
| `NHEJ_INDEL_%` | Classical NHEJ only |
| `MMEJ_INDEL_%` | MMEJ only |
| `NHEJ_MMEJ_total_%` | Combined NHEJ + MMEJ rate |
| `WT_%` | Wild-type percentage |
| `NON_DONOR_SNV_%` | Non-donor SNVs |
| `UNCLASSIFIED_%` | Unclassified reads |
| `Edited_total_%` | All edited reads (excludes WT) |

### Pooled Summary (Technical Replicates)

For experiments with technical replicates, use `trace aggregate` to generate pooled statistics:

```
bio_sample  quality_flag  n_reps  hdr_total_pct_mean  hdr_total_pct_sem  nhej_mmej_total_pct_mean
────────────────────────────────────────────────────────────────────────────────────────────────────
gene_A      good          3       45.2                1.3                12.1
gene_B      good          3       52.8                0.9                8.4
gene_C      good          3       38.1                2.1                18.7
```

The `quality_flag` column indicates:
- `good`: All replicates passed QC
- `one_low_read_removed`: One replicate had low reads and was excluded
- `all_low_reads`: All replicates had low read counts (use with caution)

## How It Works

### Triple-Aligner Consensus

TRACE uses three independent aligners to maximize accuracy:

```
                    ┌─────────────┐
     Raw Reads ────►│   BWA-MEM   │────┐
                    └─────────────┘    │
                    ┌─────────────┐    │     ┌───────────────┐     ┌────────────┐
     Raw Reads ────►│   BBMap     │────┼────►│   Consensus   │────►│  Classify  │
                    └─────────────┘    │     └───────────────┘     └────────────┘
                    ┌─────────────┐    │
     Raw Reads ────►│  minimap2   │────┘
                    └─────────────┘
```

Each aligner has different strengths—BWA-MEM for accuracy, BBMap for gapped alignments, minimap2 for speed. TRACE takes the consensus to reduce aligner-specific artifacts.

### Edit Distance Classification

TRACE classifies reads by measuring edit distance to both reference and donor sequences:

1. **Align** read to reference amplicon
2. **Extract** mismatches in core edit region (±30bp from cut site)
3. **Check** if each mismatch brings sequence closer to donor template
4. **Classify** based on SNV pattern and presence of indels

This handles complex cases where aligners represent clustered SNVs as indels.

### Automatic Detection

TRACE auto-detects library characteristics:

| Detection | Method |
|-----------|--------|
| **Library type** | TruSeq (fixed primers) vs Tn5 (tagmented) based on read start clustering |
| **UMI presence** | Sequence diversity analysis at read starts |
| **Read overlap** | Determines if R1/R2 can be merged |
| **PCR duplicates** | UMI-based (TruSeq) or position-based (Tn5) deduplication |

## Installation

### pip (Python package only)

```bash
pip install trace-crispr
```

### conda (includes aligners)

```bash
conda install -c bioconda -c conda-forge trace-crispr
```

This installs BWA, BBMap, minimap2, and samtools automatically.

### Development

```bash
git clone https://github.com/k-roy/TRACE.git
cd TRACE
pip install -e ".[dev]"
```

## Usage

### Check Locus Configuration

Before running, verify TRACE correctly detects your editing design:

```bash
trace info \
  --reference amplicon.fasta \
  --hdr-template donor.fasta \
  --guide GCTGAAGCACTGCACGCCGT
```

Output:
```
============================================================
TRACE Analysis Configuration
============================================================

Reference sequence: 231 bp
HDR template: 127 bp
  - Template aligns at position 53 in reference

Donor template analysis:
  - Left homology arm: positions 53-124 (72 bp)
  - Right homology arm: positions 128-179 (52 bp)

  Edits detected (2 total):
    * Position 125: T -> A (substitution)
    * Position 127: G -> A (substitution)

Guide analysis:
  - Guide: GCTGAAGCACTGCACGCCGT
  - Target: positions 105-124 (+ strand)
  - PAM: TGG at positions 125-127
  - Cleavage site: position 122
```

### Sample Key Format

Create a TSV file with sample information:

```
sample_id    r1_path                  r2_path                  condition
sample_01    /path/to/S1_R1.fastq.gz  /path/to/S1_R2.fastq.gz  treatment
sample_02    /path/to/S2_R1.fastq.gz  /path/to/S2_R2.fastq.gz  treatment
sample_03    /path/to/S3_R1.fastq.gz  /path/to/S3_R2.fastq.gz  control
```

### Per-Sample Sequences

For experiments with different guides/donors per sample:

```
sample_id  r1_path      r2_path      reference      guide                   hdr_template
sample_01  S1_R1.fq.gz  S1_R2.fq.gz  locus1.fasta   AGAGAAACACACTGTACTCCGT  donor1.fasta
sample_02  S2_R1.fq.gz  S2_R2.fq.gz  locus2.fasta   TTGGTTACAACTCTGACCCA    donor2.fasta
```

```bash
trace run --sample-key manifest.tsv --output results/
```

### Multi-Template Analysis (Barcode Screening)

For experiments with multiple possible HDR templates:

```bash
trace multi-template \
  --reference amplicon.fasta \
  --hdr-templates all_barcodes.fasta \
  --guide GAGTCCGAGCAGAAGAAGAA \
  --sample-key samples.tsv \
  --output results/
```

### Cas12a Support

```bash
trace run \
  --reference amplicon.fasta \
  --hdr-template donor.fasta \
  --guide GCTGAAGCACTGCACGCCGTAA \
  --nuclease cas12a \
  --sample-key samples.tsv \
  --output results/
```

| Nuclease | PAM | Cleavage |
|----------|-----|----------|
| Cas9 | NGG (3' of guide) | Blunt, 3bp from PAM |
| Cas12a | TTTN (5' of guide) | Staggered, 18-23bp from PAM |

### Snakemake Workflow

For large-scale processing:

```bash
snakemake --configfile config.yaml --cores 24
```

See [workflow/README.md](workflow/README.md) for details.

## CLI Commands Reference

TRACE provides the following commands:

| Command | Description |
|---------|-------------|
| `trace run` | Main analysis pipeline for single or multiple samples |
| `trace info` | Display locus configuration (reference, donor, guide alignment) |
| `trace classify` | Re-classify reads from an existing BAM file |
| `trace classify-batch` | Batch re-classification on multiple BAM files |
| `trace multi-template` | Analyze samples with multiple possible HDR templates |
| `trace aggregate` | Aggregate results across technical replicates |
| `trace extract-kmers` | Extract k-mers for contamination filtering |
| `trace generate-manifest` | Generate a sample manifest from a directory of FASTQs |
| `trace generate-templates` | Generate multi-template reference from barcode list |
| `trace init` | Create a template configuration file |

Use `trace <command> --help` for detailed options.

### Example: Aggregating Replicates

```bash
# After running trace on all samples
trace aggregate \
  --input results/per_sample_editing_outcomes_all_methods.tsv \
  --sample-key samples.tsv \
  --group-by condition \
  --output results/aggregated_by_condition.tsv
```

## Analysis Module

Compare editing outcomes across conditions with statistical testing:

```python
from trace_crispr.analysis import (
    compare_metric_by_condition,
    results_to_dataframe,
    plot_condition_comparison,
)

# Compare HDR rates
comparisons = compare_metric_by_condition(
    results, samples,
    condition_col='treatment',
    metric='dedup_hdr_pct',
    base_condition='control'
)

# View statistics
print(comparisons.to_dataframe())
```

```
     condition        metric  mean   std   p_value  significance
0  treatment_A  dedup_hdr_pct  25.4  2.63   0.0003          ***
1  treatment_B  dedup_hdr_pct  12.1  1.89   0.4521           ns
```

### Visualization

```python
fig = plot_condition_comparison(
    stats, comparisons,
    base_condition='control',
    title='HDR Rate by Treatment',
    ylabel='HDR Rate (%)'
)
fig.savefig('hdr_comparison.png', dpi=150)
```

Install visualization dependencies:

```bash
pip install trace-crispr[visualization]
```

## Dependencies

**Python (core):** click, pysam, pandas, numpy, pyyaml, rapidfuzz, tqdm

**Python (visualization):** matplotlib, seaborn, scipy

**External (via conda):** bwa, bbmap, minimap2, samtools

## Citation

If you use TRACE in your research, please cite:

> Roy, K.R. et al. (2026). TRACE: Triple-aligner Read Analysis for CRISPR Editing. *In preparation.*

## Author

Kevin R. Roy (kevinrjroy@gmail.com)

## License

MIT
