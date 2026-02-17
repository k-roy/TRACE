# TRACE

**T**riple-aligner **R**ead **A**nalysis for **C**RISPR **E**diting

TRACE is a comprehensive tool for quantifying CRISPR editing outcomes from amplicon sequencing data. It combines multiple alignment strategies with k-mer classification to provide robust, accurate measurements of HDR, NHEJ, and other editing outcomes.

## Features

- **Triple-aligner consensus**: Uses BWA-MEM, BBMap, and minimap2 for robust alignment
- **Flexible input**: Accepts DNA sequences directly or FASTA file paths
- **Automatic inference**: Detects PAM, cleavage site, homology arms, and edits from sequences
- **Large edit support**: Handles insertions up to 50+ bp with automatic k-mer size adjustment
- **K-mer classification**: Fast pre-alignment HDR/WT detection (auto-sizes k-mers based on edit)
- **Multi-nuclease support**: Cas9 and Cas12a (Cpf1) with correct cleavage geometry
- **Auto-detection**: Library type (TruSeq/Tn5), UMI presence, read merging need
- **PCR deduplication**: Automatic UMI-based (TruSeq) or position-based (Tn5) deduplication
- **CRISPResso2 integration**: Validation with standard CRISPR analysis tool

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

Reference sequence: 250 bp
HDR template: 150 bp
  - Template aligns at position 36 in reference

Donor template analysis:
  - Left homology arm: positions 36-110 on reference (75 bp)
  - Right homology arm: positions 113-185 on reference (73 bp)

  Edits detected (2 total):
    * Position 111: C -> G (substitution)
    * Position 112: C -> T (substitution)

Guide analysis:
  - Guide sequence: GCTGAAGCACTGCACGCCGT
  - Guide targets: positions 113-132 on reference (- strand)
  - PAM: CCC at positions 110-112 on reference
  - Cleavage site: position 116 on reference
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

### Per-sample locus sequences

For multiplexed experiments with different target loci, you can specify `reference`, `hdr_template`, and `guide` per-sample in the sample key TSV:

```
sample_id	r1_path	r2_path	condition	reference	hdr_template	guide
sample_1	/path/S1_R1.fq.gz	/path/S1_R2.fq.gz	treatment
sample_2	/path/S2_R1.fq.gz	/path/S2_R2.fq.gz	control
sample_3	/path/S3_R1.fq.gz	/path/S3_R2.fq.gz	treatment	locus2.fasta	locus2_hdr.fasta	ACGTACGTACGTACGTACGT
sample_4	/path/S4_R1.fq.gz	/path/S4_R2.fq.gz	control	locus2.fasta	locus2_hdr.fasta	ACGTACGTACGTACGTACGT
```

- **Empty values**: Use CLI defaults (`--reference`, `--hdr-template`, `--guide`)
- **Filled values**: Override CLI defaults for that sample
- Values can be DNA sequences or FASTA file paths

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

The main output is a TSV file with per-sample editing outcomes:

| Column | Description |
|--------|-------------|
| sample | Sample ID |
| classifiable_reads | Total classifiable reads |
| duplicate_rate | PCR duplicate rate |
| WT_% | Wild-type % |
| HDR_% | HDR % |
| NHEJ_% | NHEJ % |
| LgDel_% | Large deletion % |
| kmer_hdr_rate | K-mer method HDR rate |
| crispresso_hdr_rate | CRISPResso2 HDR rate |
| crispresso_indel_rate | CRISPResso2 indel rate |

For Tn5/tagmented data or TruSeq amplicons with UMIs, TRACE will report on the PCR duplication rate and automatically perform deduplication:
- **TruSeq with UMIs**: Pre-alignment UMI-based deduplication
- **Tn5**: Post-alignment position-based deduplication

## Dependencies

### Python
- click>=8.0
- pysam>=0.20
- pandas>=1.5
- numpy>=1.20
- pyyaml>=6.0
- rapidfuzz>=3.0
- tqdm>=4.60

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
