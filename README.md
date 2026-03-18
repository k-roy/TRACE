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

TRACE produces a per-sample summary table with editing outcomes:

```
sample          total_reads  HDR_%   NHEJ_%  WT_%   LgDel_%
─────────────────────────────────────────────────────────────
sample_01       125,432      45.2    12.3    38.1   4.4
sample_02       98,765       52.1    8.7     35.8   3.4
sample_03       112,890      48.9    15.2    31.2   4.7
```

**Classification categories:**
| Outcome | Description |
|---------|-------------|
| `HDR_COMPLETE` | All donor SNVs integrated, no additional edits |
| `HDR_PARTIAL` | Subset of donor SNVs, no additional edits |
| `NHEJ` | Indels near cut site, no donor sequence |
| `MIXED` | Both donor SNVs and non-donor edits |
| `WT` | No edits detected |
| `LARGE_DELETION` | Deletions spanning >50bp |

### Conversion Tract Analysis

For HDR samples, TRACE outputs per-SNV integration frequencies showing how donor sequence propagates from the cut site:

```
position  distance_to_cut  ref  donor  frequency  count
────────────────────────────────────────────────────────
125       -3               T    A      0.95       1,234
127       -1               G    C      0.92       1,198
130       +2               A    T      0.87       1,132
135       +7               C    G      0.71       923
142       +14              T    A      0.45       585
```

**Interpretation:** SNVs near the cut site show high integration (>90%), decreasing with distance. This reveals conversion tract length and strand bias in HDR repair.

### Pooled Summary (Technical Replicates)

For experiments with technical replicates, TRACE generates a pooled summary with weighted statistics:

```
bio_sample  n_reps  HDR_pct_mean  HDR_pct_sem  NHEJ_pct_mean  total_reads_sum
──────────────────────────────────────────────────────────────────────────────
gene_A      3       45.2          1.3          12.1           356,087
gene_B      3       52.8          0.9          8.4            312,456
gene_C      3       38.1          2.1          18.7           298,234
```

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
