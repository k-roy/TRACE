# TRACE Snakemake Workflow

This directory contains a Snakemake workflow for running TRACE analysis with automatic dependency tracking and resumability.

## Features

- **Automatic skip detection**: Only reprocesses samples when inputs change
- **Parallel execution**: Processes multiple samples simultaneously
- **Per-sample sequences**: Supports different guide/donor/reference per sample
- **Resumable**: Can restart from any point if interrupted
- **Modular**: Separate rules for trim, align, classify

## Installation

Install Snakemake:
```bash
conda install -c bioconda snakemake
# or
pip install snakemake
```

## Quick Start

### 1. Prepare your manifest

Create a TSV file with per-sample sequences (optional columns):

```tsv
sample_id	r1_path	r2_path	reference	guide	hdr_template
sample_1	S1_R1.fq.gz	S1_R2.fq.gz	ATCG...	AGAGAAACACACTGTACTCCGT	ATCG...
sample_2	S2_R1.fq.gz	S2_R2.fq.gz	ref.fasta	TTGGTTACAACTCTGACCCA	hdr.fasta
```

Or use global defaults in `config.yaml`.

### 2. Configure workflow

Edit `config.yaml`:
```yaml
manifest: "samples.tsv"
output_dir: "results/trace"
threads_per_sample: 4
default_reference: "reference.fasta"  # If not in manifest
default_guide: "GCTGAAGCACTGCACGCCGT"
default_hdr_template: "hdr_template.fasta"
```

### 3. Run workflow

**Dry run** (see what will be executed):
```bash
snakemake --cores 24 --configfile config.yaml -n
```

**Full run**:
```bash
snakemake --cores 24 --configfile config.yaml
```

**Resume from checkpoint**:
```bash
snakemake --cores 24 --configfile config.yaml
```
Snakemake automatically detects existing outputs and skips completed steps.

## Workflow Steps

```
Raw FASTQ (R1, R2)
    ↓
trim_adapters → Trimmed FASTQ
    ↓
align_bwa + align_bbmap + align_minimap2 → Triple alignment (BWA, BBMap, minimap2)
    ↓
classify_reads (uses per-sample guide/donor) → classification.tsv
    ↓
aggregate_results → all_samples_results.tsv
```

**Triple-aligner consensus**: Following TRACE's core design, the workflow runs all three aligners (BWA-MEM, BBMap, minimap2) and uses the first successful aligner for classification. This provides robust alignment across diverse read types and reference sequences.

## Advanced Usage

### Only classify (reuse existing alignments)

If you already have trimmed FASTQs and alignments:
```bash
snakemake --cores 24 --forcerun classify_reads aggregate_results
```

### Run specific samples
```bash
snakemake results/trace/sample_1/classification.tsv --cores 4
```

### Migrate existing TRACE outputs

If you have existing TRACE outputs and want to reuse alignments with the Snakemake workflow:

```python
#!/usr/bin/env python3
"""Migrate existing TRACE outputs to Snakemake directory structure"""
from pathlib import Path
import pandas as pd

manifest = pd.read_csv("trace_manifest.tsv", sep='\t')
old_dir = Path("results/trace_original")
new_dir = Path("results/trace_snakemake")

for _, row in manifest.iterrows():
    sample_id = row['sample_id']

    # Symlink trimmed FASTQs
    old_r1 = old_dir / sample_id / "trimmed" / "R1_trimmed.fastq.gz"
    old_r2 = old_dir / sample_id / "trimmed" / "R2_trimmed.fastq.gz"
    new_r1 = new_dir / sample_id / "trimmed" / "R1_trimmed.fastq.gz"
    new_r2 = new_dir / sample_id / "trimmed" / "R2_trimmed.fastq.gz"

    if old_r1.exists():
        new_r1.parent.mkdir(parents=True, exist_ok=True)
        new_r1.symlink_to(old_r1.resolve())
        new_r2.symlink_to(old_r2.resolve())

    # Symlink BAM files (all three aligners)
    for aligner in ['bwa', 'bbmap', 'minimap2']:
        old_bam = old_dir / sample_id / "alignments" / f"{aligner}.bam"
        new_bam = new_dir / sample_id / "alignments" / f"{aligner}.bam"

        if old_bam.exists():
            new_bam.parent.mkdir(parents=True, exist_ok=True)
            new_bam.symlink_to(old_bam.resolve())

print(f"Migrated outputs from {old_dir} to {new_dir}")
print("Now run: snakemake --forcerun classify_reads aggregate_results --cores 24")
```

Then force reclassification with updated sequences:
```bash
snakemake --forcerun classify_reads aggregate_results --cores 24
```

### Visualize workflow
```bash
snakemake --dag | dot -Tpng > workflow.png
```

### Clean outputs
```bash
snakemake --delete-all-output
```

## Output Structure

```
results/trace/
├── reference.fasta                    # Reference for alignment
├── sample_1/
│   ├── trimmed/
│   │   ├── R1_trimmed.fastq.gz
│   │   └── R2_trimmed.fastq.gz
│   ├── alignments/
│   │   ├── bwa.bam                    # BWA-MEM alignment
│   │   ├── bbmap.bam                  # BBMap alignment
│   │   └── minimap2.bam               # minimap2 alignment
│   ├── classification.tsv             # Per-read classifications
│   └── logs/
│       ├── trim.log
│       ├── bwa.log
│       ├── bbmap.log
│       ├── minimap2.log
│       └── classify.log
├── sample_2/
│   └── ...
└── summary/
    └── all_samples_results.tsv        # Aggregated results
```

**Classification output columns:**
- `read_name`: Read identifier
- `outcome`: HDR, NHEJ, WT, or LARGE_DELETION
- `hdr_fraction`: Fraction of HDR signature positions matched (0.0-1.0)

**Summary output columns:**
- `sample_id`: Sample identifier
- `total_reads`: Number of classified reads
- `HDR`, `NHEJ`, `WT`, `LARGE_DELETION`: Read counts by outcome
- `HDR_pct`: HDR percentage

## Integration with Per-Sample Sequences (TRACE v0.5.0)

The workflow automatically uses per-sample `reference`, `guide`, and `hdr_template` from the manifest:

- If columns exist in manifest → use per-sample values
- If columns are empty → use config defaults
- Mix and match: some samples with custom sequences, others with defaults

**Key advantage**: When you update guide/donor sequences in the manifest, Snakemake only reruns classification (not expensive trimming/alignment steps).

**Implementation note**: The workflow uses Snakemake's `script:` directive with external Python files ([scripts/classify_sample.py](scripts/classify_sample.py), [scripts/aggregate_results.py](scripts/aggregate_results.py)) to avoid subprocess import issues. This ensures TRACE modules are correctly imported in Snakemake's isolated rule execution environment.

## Troubleshooting

### Workflow won't resume

Force rerun of specific rule:
```bash
snakemake --forcerun <rule_name> --cores 24
```

### Need to change alignment parameters

Edit the rule in `Snakefile`, then:
```bash
snakemake --forcerun align_bwa --cores 24
```

### Check why a rule will run
```bash
snakemake <target> --reason
```

## Performance Tips

1. **Set threads wisely**: `threads_per_sample × concurrent_samples ≤ total_cores`
   ```bash
   # Example: 24 cores, 4 threads per sample → max 6 concurrent samples
   snakemake --cores 24  # Snakemake auto-schedules
   ```

2. **Use local temp for I/O-heavy steps**:
   ```python
   temp(OUTPUT_DIR / "{sample}" / "trimmed" / "R1.fastq.gz")
   ```

3. **Profile to find bottlenecks**:
   ```bash
   snakemake --cores 24 --profile profile/
   ```
