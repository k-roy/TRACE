# TRACE

**T**riple-aligner **R**ead **A**nalysis for **C**RISPR **E**diting

## Features

- **Triple-aligner consensus**: Uses BWA-MEM, BBMap, and minimap2 for robust alignment
- **Automatic inference**: Detects PAM, cleavage site, homology arms, and edits from sequences
- **K-mer classification**: Fast pre-alignment HDR/WT detection using 12-mers
- **Multi-nuclease support**: Cas9 and Cas12a (Cpf1) with correct cleavage geometry
- **Auto-detection**: Library type (TruSeq/Tn5), read merging need, CRISPResso mode
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
git clone https://github.com/k-roy/trace.git
cd trace
pip install -e ".[dev]"
```

## Quick Start

### Minimal run (3 required inputs)

```bash
trace run \
  --reference amplicon.fasta \
  --hdr-template hdr_template.fasta \
  --guide GCTGAAGCACTGCACGCCGT \
  --r1 sample_R1.fastq.gz \
  --r2 sample_R2.fastq.gz \
  --output results/
```

### Check locus configuration without running

```bash
trace info \
  --reference amplicon.fasta \
  --hdr-template hdr_template.fasta \
  --guide GCTGAAGCACTGCACGCCGT
```

This will print:

```
=== TRACE Analysis Configuration ===

Reference sequence: 500 bp
HDR template: 500 bp

Donor template analysis:
  - Left homology arm: positions 1-245 on reference (245 bp)
  - Right homology arm: positions 255-500 on reference (245 bp)
  - Donor edits detected at positions: 246, 247 on reference
    * Position 246: C → G (PAM-silencing mutation)
    * Position 247: C → T (chromophore Y66H mutation)

Guide analysis:
  - Guide sequence: GCTGAAGCACTGCACGCCGT
  - Guide targets: positions 248-267 on reference (- strand)
  - PAM: GGG at positions 245-247 on reference
  - Cleavage site: position 248 on reference
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
| duplicate_rate | PCR duplicate rate (Tn5) |
| Dedup_WT_% | Wild-type % (deduplicated) |
| Dedup_HDR_% | HDR % (deduplicated) |
| Dedup_NHEJ_% | NHEJ % (deduplicated) |
| Dedup_LgDel_% | Large deletion % |
| kmer_hdr_rate | K-mer method HDR rate |
| crispresso_hdr_rate | CRISPResso2 HDR rate |

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
