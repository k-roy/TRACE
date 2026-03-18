# TRACE

**T**riple-aligner **R**ead **A**nalysis for **C**RISPR **E**diting

TRACE is a comprehensive tool for quantifying CRISPR editing outcomes from amplicon sequencing data. It combines multiple alignment strategies with k-mer classification to provide robust, accurate measurements of HDR, NHEJ, and other editing outcomes.

## Features

- **Triple-aligner consensus**: Uses BWA-MEM, BBMap, and minimap2 for robust alignment
- **Flexible input**: Accepts DNA sequences directly or FASTA file paths
- **Automatic inference**: Detects PAM, cleavage site, homology arms, and edits from sequences
- **Large edit support**: Handles insertions up to 50+ bp with automatic k-mer size adjustment
- **K-mer classification**: Fast pre-alignment HDR/WT detection (auto-sizes k-mers based on edit)
- **Barcode-optimized k-mers**: Auto-detects barcode-style templates and generates discriminating k-mers
- **Multi-nuclease support**: Cas9 and Cas12a (Cpf1) with correct cleavage geometry
- **Robust UMI detection**: 3-pass algorithm handles variable primer quality and low-signal libraries
- **Auto-detection**: Library type (TruSeq/Tn5), UMI presence, read merging need
- **PCR deduplication**: Automatic UMI-based (TruSeq) or position-based (Tn5) deduplication
- **CRISPResso2 integration**: Validation with standard CRISPR analysis tool

## Recent Updates

### Version 0.6.0 (2026-03-07)

**New Features:**
- **Edit distance HDR detection**: New default classification method that handles clustered SNVs correctly
  - Checks if each read mismatch brings the sequence closer to the donor template
  - Filters to "core edit region" (±30bp from cut site) to exclude flanking sequence variations
  - Tracks per-SNV integration for conversion tract analysis
  - Handles cases where aligners represent clustered SNVs as indels

- **Conversion tract analysis**: Per-SNV integration frequencies by distance from cut site
  - Outputs `hdr_snv_detail.tsv` with integration frequency at each donor SNV position
  - Enables analysis of how far HDR extends from the cut site
  - Useful for understanding donor template design and nickase strand effects

**Classification outcomes:**
- `HDR_COMPLETE`: All expected donor SNVs present, no additional edits
- `HDR_PARTIAL`: Subset of donor SNVs present, no additional edits
- `NHEJ`: Indels near cut site, no donor SNVs
- `MIXED`: Both donor SNVs and non-donor edits (potential NHEJ + partial HDR)
- `WT`: No edits detected

### Version 0.5.0 (2026-03-05)

**New Features:**
- **Per-sample guide/donor/reference sequences**: Each sample can now use different guide RNAs, donor templates, and reference sequences from the manifest file
- **Snakemake workflow**: Official workflow for parallel processing with automatic dependency tracking and resumability
- **Flexible sequence input**: Specify per-sample sequences in manifest columns or use CLI defaults for all samples
- **Locus config caching**: Per-sample configurations are cached by sequence content to optimize processing of biological replicates

**Usage:**
```bash
# All samples with per-sample sequences (no CLI defaults needed)
trace run --sample-key manifest.tsv --output results/

# Mix of per-sample and default sequences
trace run \
  --reference "ATCG..." \
  --guide "AGAG..." \
  --hdr-template "ATCG..." \
  --sample-key manifest.tsv \
  --output results/
```

**Snakemake Workflow:**
```bash
# Full pipeline (trim → triple-align → classify)
snakemake --configfile config.yaml --cores 24

# Classification only (if alignments exist)
snakemake --forcerun classify_reads aggregate_results --cores 24
```

The Snakemake workflow runs all three aligners (BWA-MEM, BBMap, minimap2) and uses the first successful aligner for classification, matching the core TRACE triple-aligner approach. See [workflow/README.md](workflow/README.md) for detailed documentation.

### Version 0.4.0 (2026-02-20)

**New Features:**
- **Alignment-only classification (now default)**: Pure alignment-based classification using multi-reference FASTA (WT + all HDR variants). More accurate than k-mer classification, especially for experiments with many similar barcodes.
- **Full-length HDR builder**: Automatically builds full amplicon sequences from short donor templates using homology arm detection
- **Global sequence deduplication**: Filter low-count sequences across all samples with `min_global_count` parameter
- **Primary-only alignments**: When aligning to multi-reference FASTA, secondary alignments are suppressed by default

**Breaking Changes:**
- `alignment_only` parameter now defaults to `True`. Set `alignment_only=False` for legacy k-mer mode.

**Bug Fixes:**
- Fixed BBMap subprocess deadlock caused by stderr buffer overflow
- Fixed version mismatch between `__init__.py` and `pyproject.toml`

See [CHANGELOG.md](CHANGELOG.md) for full details.

### Version 0.3.1 (2026-02-17)

**Critical Bug Fixes:**
- **3-pass UMI detection algorithm**: Dramatically improves merge rates for libraries with weak primer signals
  - Pass 1: Strong signal detection (>50% consensus) - high confidence
  - Pass 2: Weak signal detection (>30% consensus, ≥4bp UMI) - medium confidence
  - Pass 3: Jump detection with 6bp fallback - handles poor quality data
  - **Impact**: Merge rate improved from 35% to 92.5% on HEK293 EMX1 test dataset (80 samples)

**Optimizations:**
- **Barcode-style template auto-detection**: Automatically identifies when templates share homology arms with different barcodes
- **Optimized k-mer generation**: Generates k-mers that span barcode boundaries for better discrimination
- **Multi-template batch processing**: ~100x faster classification through global sequence deduplication

**Dependencies:**
- Added CRISPResso2 to optional validation dependencies (`pip install trace-crispr[validation]`)

## Installation

### pip (Python package only)

```bash
pip install trace-crispr
```

### conda (includes external aligners)

```bash
conda install -c bioconda -c conda-forge trace-crispr
```

### Development installation

```bash
git clone https://github.com/k-roy/TRACE.git
cd TRACE
pip install -e ".[dev]"
```

## Quick Start

TRACE accepts sequences as either **DNA strings** or **FASTA file paths**.

### Example 1: Using FASTA files

```bash
trace run \
  --reference amplicon.fasta \
  --hdr-template hdr_template.fasta \
  --guide GCTGAAGCACTGCACGCCGT \
  --r1 sample_R1.fastq.gz \
  --r2 sample_R2.fastq.gz \
  --output results/
```

### Example 2: Using DNA sequences directly

The reference amplicon will typically be longer than the HDR template so that the primers specifically amplify the target locus. This example shows a 250 bp reference sequence where the donor template is 150 bp long with the designed edit in the middle.

```bash
# Reference amplicon (250 bp) - includes flanking regions
# Guide sequence shown in lowercase for illustration
REF="ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\
gctgaagcactgcacgccgttgg\
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"

# HDR template (150 bp) - centered on edit site
# Guide in lowercase, designed edit (T->A) shown in UPPERCASE
HDR="ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\
gctgaagcactgcacgccgtAga\
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"

trace run \
  -r "$REF" \
  -h "$HDR" \
  -g GCTGAAGCACTGCACGCCGT \
  --r1 sample_R1.fastq.gz \
  --output results/
```

**Note:** The guide and edit are shown here in lowercase/uppercase for illustrative purposes. This formatting is not necessary - TRACE will automatically detect guide and edit positions from the sequences.

### Check locus configuration without running

```bash
trace info \
  --reference amplicon.fasta \
  --hdr-template hdr_template.fasta \
  --guide GCTGAAGCACTGCACGCCGT
```

This will print:

```
============================================================
=== TRACE Analysis Configuration ===
============================================================

Reference sequence: 231 bp
HDR template: 127 bp
  - Template aligns at position 53 in reference

Donor template analysis:
  - Left homology arm: positions 53-124 on reference (72 bp)
  - Right homology arm: positions 128-179 on reference (52 bp)

  Edits detected (2 total):
    * Position 125: T -> A (substitution)
    * Position 127: G -> A (substitution)

Guide analysis:
  - Guide sequence: GCTGAAGCACTGCACGCCGT
  - Guide targets: positions 105-124 on reference (+ strand)
  - PAM: TGG at positions 125-127 on reference
  - Cleavage site: position 122 on reference
```

### Multiple samples

Create a sample key TSV:

```
sample_id	r1_path	r2_path	condition
sample_1	/path/to/S1_R1.fastq.gz	/path/to/S1_R2.fastq.gz	treatment
sample_2	/path/to/S2_R1.fastq.gz	/path/to/S2_R2.fastq.gz	control
```

Then run:

```bash
trace run \
  --reference amplicon.fasta \
  --hdr-template hdr_template.fasta \
  --guide GCTGAAGCACTGCACGCCGT \
  --sample-key samples.tsv \
  --output results/ \
  --threads 16
```

### Per-sample locus sequences (NEW in v0.5.0)

**TRACE now supports heterogeneous experiments where different samples use different guide RNAs, donor templates, or reference sequences.**

For experiments with different target loci, editing strategies, or guides across samples, you can specify `reference`, `hdr_template`, `guide`, and optionally `nuclease` per-sample in the sample key TSV:

```
sample_id	r1_path	r2_path	reference	guide	hdr_template
sample_1	S1_R1.fq.gz	S1_R2.fq.gz	ATCG...	AGAGAAACACACTGTACTCCGT	ATCG...
sample_2	S2_R1.fq.gz	S2_R2.fq.gz	ATCG...	TTGGTTACAACTCTGACCCA	ATCG...
sample_3	S3_R1.fq.gz	S3_R2.fq.gz	locus2.fasta	ACGTACGTACGTACGTACGT	locus2_hdr.fasta
```

**Flexible usage patterns:**

#### Option 1: All samples with per-sample sequences
```bash
trace run \
  --sample-key manifest_with_sequences.tsv \
  --output results/ \
  --threads 16
```
No CLI defaults needed if all samples have complete sequence columns.

#### Option 2: Mix of per-sample and default sequences
```bash
trace run \
  --reference default_ref.fasta \      # Default for samples without custom reference
  --guide GCTGAAGCACTGCACGCCGT \       # Default for samples without custom guide
  --hdr-template default_hdr.fasta \   # Default for samples without custom template
  --sample-key manifest.tsv \
  --output results/
```
Samples with per-sample columns use those sequences; others use CLI defaults.

**Sequence input formats:**
- **DNA strings**: `ATCGATCGATCG...` (case-insensitive, converted to uppercase)
- **FASTA files**: Path to `.fa` or `.fasta` file (first sequence used)

**Performance note:**
TRACE caches per-sample locus configurations. Samples with identical sequences (e.g., biological replicates) share the same locus config, k-mer classifier, and HDR signature for efficiency.

**Example use case:**
Testing sgRNA_1 vs sgRNA_2 on the same locus, or comparing 125bp vs 346bp donors across samples in a single analysis run.

### Multi-Template Analysis (Barcode Screening)

For barcode screening experiments where multiple HDR templates (barcodes) are possible, use the `multi-template` command:

```bash
# Generate HDR templates FASTA from keyfiles
trace generate-templates \
  --sample-key keyfiles/sample_key.tsv \
  --seq-ref keyfiles/guide_donor_and_reference_info.tsv \
  --output templates/hdr_templates.fasta

# Generate sample manifest from keyfiles
trace generate-manifest \
  --sample-key keyfiles/sample_key.tsv \
  --plate-key keyfiles/plate_key.tsv \
  --raw-data-dir raw_data/ \
  --output trace_sample_key.tsv

# Run multi-template analysis
trace multi-template \
  --reference templates/reference.fasta \
  --hdr-templates templates/hdr_templates.fasta \
  --guide GAGTCCGAGCAGAAGAAGAA \
  --sample-key trace_sample_key.tsv \
  --output results/ \
  --threads 16
```

**Output tables:**
- `per_sample_editing_outcomes_all_methods.tsv` - Summary per sample
- `per_sample_per_template_outcomes.tsv` - Granular per-template results

**Features:**
- Detects which barcode/template is present in each read
- Purity checking (detects unexpected barcodes)
- AMBIGUOUS category for reads matching multiple barcodes
- Expected template validation via `expected_barcode` column
- Parallel processing with configurable threads

### Using Cas12a

```bash
trace run \
  --reference amplicon.fasta \
  --hdr-template hdr_template.fasta \
  --guide GCTGAAGCACTGCACGCCGTAA \
  --nuclease cas12a \
  --sample-key samples.tsv \
  --output results/
```

## Auto-Detection

TRACE automatically detects library characteristics to optimize analysis:

### Library Type Detection

TRACE distinguishes between **TruSeq** (fixed-target amplicon) and **Tn5** (tagmented locus) libraries by analyzing read alignment positions:

- **TruSeq**: Reads cluster at fixed start positions (primer binding sites)
- **Tn5**: Reads have scattered start positions (random Tn5 cutting)

```
Auto-detection results:
  - Library type: TruSeq (100% of reads cluster at fixed start position)
  - UMI detection: UMIs of length 6 bp detected
    --> Entering PCR deduplication mode...
```

### UMI Detection

For TruSeq libraries, TRACE detects UMIs (Unique Molecular Identifiers) by analyzing sequence diversity at read starts:

- High diversity region = UMI
- Low diversity region = primer sequence
- Automatically determines UMI length (typically 4-12 bp)

### Preprocessing Modes

Based on detection results, TRACE automatically selects the optimal preprocessing workflow:

| Library | UMIs | Overlap (>=15bp) | Preprocessing | Output |
|---------|------|------------------|---------------|--------|
| TruSeq | Yes | Yes | dedup → trim → merge → collapse | merged FASTQ |
| TruSeq | Yes | No | dedup → trim | paired FASTQs |
| TruSeq | No | Yes | trim → merge → collapse | merged FASTQ |
| TruSeq | No | No | trim | paired FASTQs |
| Tn5 | N/A | Yes | trim → merge → collapse → align → position dedup | merged FASTQ |
| Tn5 | N/A | No | trim → align → position dedup | paired FASTQs |

Example output:
```
Auto-detection results:
  - Library type: TruSeq (100% of reads cluster at fixed start position)
  - UMI detection: UMIs of length 6 bp detected
  - Read overlap: Enabled (~50bp overlap (25% of amplicon))
  - Preprocessing: dedup-trim-merge-collapse -> merged
    (UMI dedup -> trim -> merge -> collapse)
  - CRISPResso mode: merged (merged reads from preprocessing)
```

## Designed Edit Detection

TRACE automatically detects the edits encoded in the donor by first aligning the HDR template to the reference. TRACE then classifies the intended edit as a single-nucleotide variant (SNV), multi-nucleotide variant (MNV), insertion, or deletion.

K-mers are selected that span the designed edit and are unique to the reference and donor. For large edits, TRACE automatically increases the k-mer size to ensure reliable classification. TRACE can handle MNVs or insertions with lengths up to the read length - 50 bp (e.g., 100 bp insertions can be detected with 150 bp reads).

Example output for a 20 bp insertion:
```
Edits detected (1 total):
  * Position 125: +ATCGATCGATCGATCGATCG (20 bp insertion)

Maximum edit size: 20 bp
Recommended k-mer size: 30 bp
```

## Nuclease Support

### Cas9 (SpCas9)
- PAM: NGG (3' of protospacer)
- Cleavage: 3 bp upstream of PAM (blunt ends)

### Cas12a (LbCpf1)
- PAM: TTTN (5' of protospacer)
- Cleavage: 18-19 bp downstream on target strand, 23 bp on non-target
- Creates 4-5 nt 5' overhang (staggered cut)

## Output

### Per-sample classification

The main output is a TSV file with per-sample editing outcomes:

| Column | Description |
|--------|-------------|
| sample | Sample ID |
| classifiable_reads | Total classifiable reads |
| duplicate_rate | PCR duplicate rate |
| WT_% | Wild-type % |
| HDR_% | HDR % (complete + partial) |
| HDR_COMPLETE_% | HDR with all donor SNVs |
| HDR_PARTIAL_% | HDR with subset of donor SNVs |
| NHEJ_% | NHEJ % |
| MIXED_% | Mixed HDR + NHEJ events |
| LgDel_% | Large deletion % |

### Per-SNV integration detail (conversion tract analysis)

For samples with HDR templates, TRACE outputs `hdr_snv_detail.tsv` tracking integration at each donor SNV position:

| Column | Description |
|--------|-------------|
| position | Genomic position of SNV |
| distance_to_cut | Signed distance from cut site (negative = 5' of cut) |
| ref_base | Base in reference sequence |
| donor_base | Base in donor template |
| frequency | Integration frequency among HDR reads |
| count | Number of reads with this SNV |

**Interpretation**: This file reveals **conversion tract patterns** - how far donor sequence gets incorporated from the cut site. Typical findings:
- SNVs near the cut site have higher integration frequency
- Integration frequency decreases with distance from cut
- Nickase strand preference affects the gradient (5' vs 3' of cut)

**Note**: The frequencies are among reads already classified as HDR (complete, partial, or mixed) - not overall HDR rates.

For Tn5/tagmented data or TruSeq amplicons with UMIs, TRACE will report on the PCR duplication rate and automatically perform deduplication:
- **TruSeq with UMIs**: Pre-alignment UMI-based deduplication
- **Tn5**: Post-alignment position-based deduplication

## Analysis and Visualization

TRACE includes an analysis module for comparing editing outcomes across conditions with statistical testing and publication-quality visualizations.

### Installation

The analysis module requires additional dependencies for visualization:

```bash
pip install trace-crispr[visualization]
```

Or install scipy, matplotlib, and seaborn separately:

```bash
pip install scipy matplotlib seaborn
```

### Basic Usage

#### Compare conditions from pipeline results

```python
from trace_crispr.analysis import (
    compare_metric_by_condition,
    results_to_dataframe,
    get_condition_stats,
    plot_condition_comparison,
)

# After running the TRACE pipeline
results = pipeline.run_all(samples)

# Compare HDR rates across conditions
comparisons = compare_metric_by_condition(
    results, samples,
    condition_col='treatment',      # Column in sample metadata
    metric='dedup_hdr_pct',         # Metric to compare
    base_condition='control'        # Reference condition for t-tests
)

# View results as a DataFrame
print(comparisons.to_dataframe())
```

Output:
```
     condition base_condition        metric  condition_mean  condition_std  ...  p_value  p_adjusted significance
0  treatment_A        control  dedup_hdr_pct           25.41           2.63  ...   0.0003      0.0006          ***
1  treatment_B        control  dedup_hdr_pct           12.15           1.89  ...   0.4521      0.4521           ns
```

#### Create bar plots with replicate points

```python
# Convert results to DataFrame and get stats
df = results_to_dataframe(results, samples)
stats = get_condition_stats(df, 'treatment', 'dedup_hdr_pct')

# Create bar plot with individual points and significance stars
fig = plot_condition_comparison(
    stats, comparisons,
    base_condition='control',
    title='HDR Rate by Treatment',
    ylabel='HDR Rate (%)'
)
fig.savefig('hdr_comparison.png', dpi=150, bbox_inches='tight')
```

This creates a bar chart showing:
- Mean values as bars
- Individual replicate values as overlaid points (with jitter)
- Standard error of mean (SEM) as error bars
- Significance stars above significantly different conditions

#### Work directly with DataFrames

If you already have a DataFrame (e.g., from a previous analysis):

```python
from trace_crispr.analysis import (
    compare_dataframe_by_condition,
    get_condition_stats,
    plot_condition_comparison,
)
import pandas as pd

# Load existing data
df = pd.read_csv('editing_outcomes.tsv', sep='\t')

# Compare conditions
comparisons = compare_dataframe_by_condition(
    df,
    condition_col='treatment',
    metric='dedup_hdr_pct',
    base_condition='control'
)

# Get stats and plot
stats = get_condition_stats(df, 'treatment', 'dedup_hdr_pct')
fig = plot_condition_comparison(stats, comparisons)
```

#### Get summary statistics

```python
from trace_crispr.analysis import get_condition_summary

# Generate summary table for all metrics
summary = get_condition_summary(results, samples, condition_col='treatment')
print(summary[['condition', 'n', 'dedup_hdr_pct_mean', 'dedup_hdr_pct_sem']])
```

Output:
```
     condition  n  dedup_hdr_pct_mean  dedup_hdr_pct_sem
0      control  4               10.25               0.68
1  treatment_A  4               25.41               1.32
2  treatment_B  4               12.15               0.94
```

### Statistical Methods

- **T-test**: Welch's t-test (unequal variances) comparing each condition to the base
- **FDR correction**: Benjamini-Hochberg correction applied by default (disable with `fdr_correction=False`)
- **Significance thresholds**: `*` (p < 0.05), `**` (p < 0.01), `***` (p < 0.001)

### Available Functions

| Function | Description |
|----------|-------------|
| `compare_metric_by_condition()` | Main entry point - compare a metric across conditions from SampleResults |
| `compare_dataframe_by_condition()` | Compare conditions from an existing DataFrame |
| `get_condition_summary()` | Generate summary statistics table |
| `results_to_dataframe()` | Convert SampleResult list to DataFrame |
| `get_condition_stats()` | Calculate mean, std, sem, n for each condition |
| `compare_conditions()` | Perform statistical comparisons between conditions |
| `plot_condition_comparison()` | Bar plot with points, error bars, and significance stars |
| `plot_comparison_summary()` | Forest/bar plot of fold changes |
| `plot_replicate_correlation()` | Scatter plot comparing two metrics |
| `plot_multi_metric_comparison()` | Multi-panel comparison across metrics |

### Plot Customization

```python
fig = plot_condition_comparison(
    stats, comparisons,
    base_condition='control',
    title='HDR Rate by Treatment',
    ylabel='HDR Rate (%)',
    figsize=(12, 6),              # Figure size
    bar_color='#22c55e',          # Color for treatment bars (green)
    base_color='#888888',         # Color for base condition (gray)
    point_alpha=0.6,              # Transparency for data points
    jitter=0.15,                  # Horizontal spread for points
    show_significance=True,       # Show significance stars
    condition_order=['control', 'treatment_A', 'treatment_B'],  # Custom order
)
```

## Dependencies

### Python (Core)
- click>=8.0
- pysam>=0.20
- pandas>=1.5
- numpy>=1.20
- pyyaml>=6.0
- rapidfuzz>=3.0
- tqdm>=4.60

### Python (Visualization - optional)
- matplotlib>=3.5
- seaborn>=0.12
- scipy>=1.9

Install with: `pip install trace-crispr[visualization]`

### External tools (via conda)
- bwa>=0.7
- bbmap>=39
- minimap2>=2.24
- samtools>=1.16
- crispresso2 (optional, but enabled by default)

## Author

Kevin R. Roy

## License

MIT
