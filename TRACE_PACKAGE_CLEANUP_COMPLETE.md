# TRACE Package Version Suffix Removal - Complete

**Date:** 2026-03-10
**Status:** ✅ COMPLETE

---

## Summary

Comprehensive cleanup of the TRACE software package to remove all version suffixes (_v2, _v3) from workflow files. The package now contains only the single correct version of each file, with obsolete versions archived for reference.

---

## Changes Made

### 1. ✅ Workflow Scripts Renamed

**Directory:** `/oak/stanford/groups/larsms/Users/kevinroy/software/trace/workflow/scripts/`

| Old Name | New Name | Status |
|----------|----------|--------|
| `collapse_reads_v3.py` | `collapse_reads.py` | ✅ Renamed (current version) |
| `expand_to_samples_v2.py` | `expand_to_samples.py` | ✅ Renamed (current version) |

**Archived (obsolete versions):**
- `collapse_reads.py` (old) → `archive/obsolete_workflow_scripts/`
- `collapse_reads_v2.py` (intermediate) → `archive/obsolete_workflow_scripts/`
- `expand_to_samples.py` (old) → `archive/obsolete_workflow_scripts/`
- `expand_to_samples_v2.py.backup` (backup) → `archive/obsolete_workflow_scripts/`
- `classify_sample_v2.py` (unused) → `archive/obsolete_workflow_scripts/`

---

### 2. ✅ Snakefile Renamed

**Directory:** `/oak/stanford/groups/larsms/Users/kevinroy/software/trace/workflow/`

| Old Name | New Name | Status |
|----------|----------|--------|
| `Snakefile_cached_v2` | `Snakefile_cached` | ✅ Renamed (current version) |

**Archived:**
- `Snakefile_cached` (old) → `archive/obsolete_snakefiles/`

---

### 3. ✅ Snakefile Updated

**File:** `workflow/Snakefile_cached`

**Updated script references:**
```python
# Line 157: collapse_reads_v3.py → collapse_reads.py
script:
    str(WORKFLOW_DIR / "scripts" / "collapse_reads.py")

# Line 302: expand_to_samples_v2.py → expand_to_samples.py
script:
    str(WORKFLOW_DIR / "scripts" / "expand_to_samples.py")
```

---

### 4. ✅ Project References Updated

**HEK293 LGC Reporter Project:**

Updated 8 files to reference `Snakefile_cached` instead of `Snakefile_cached_v2`:

1. **scripts/submit_trace_workflow.sh** - Main SLURM submission script
2. **regenerate_with_fixes.sh** - Regeneration script
3. **submit_trace_full.sh** - Full dataset submission
4. **README.md** - Project documentation
5. **TRACE_WORKFLOW.md** - Workflow documentation
6. **TRACE_BUG_FIXES_SUMMARY.md** - Bug fix summary
7. **CLEANUP_COMPLETE.md** - Cleanup documentation
8. **CLEANUP_PLAN.md** - Cleanup plan

**Example change:**
```bash
# OLD:
snakemake \
    --snakefile /oak/.../trace/workflow/Snakefile_cached_v2 \
    --config manifest=keyfiles/trace_manifest_full.tsv

# NEW:
snakemake \
    --snakefile /oak/.../trace/workflow/Snakefile_cached \
    --config manifest=keyfiles/trace_manifest_full.tsv
```

---

## Current Clean Structure

```
/oak/stanford/groups/larsms/Users/kevinroy/software/trace/
├── trace_crispr/                          # Core TRACE package
│   ├── core/
│   │   ├── edit_distance_hdr.py          # HDR detection (bug fixes applied)
│   │   ├── alignment_classifier.py       # Cut-site-aware NHEJ (bug fixes applied)
│   │   └── classification.py             # Outcome classification (updated)
│   └── workflow/
│       └── (workflow utilities)
│
├── workflow/                              # Snakemake workflows
│   ├── Snakefile                         # Original per-sample workflow
│   ├── Snakefile_cached                  # ✨ Two-level cached workflow (current)
│   │
│   └── scripts/                          # ✨ All scripts renamed (no version suffixes)
│       ├── collapse_reads.py             # Current version (was v3)
│       ├── expand_to_samples.py          # Current version (was v2)
│       ├── classify_per_guide_donor.py   # Level 2 classification
│       ├── extract_global_unique_seqs.py # Level 1 caching
│       ├── align_unique_sequences.py     # Alignment
│       ├── extract_alignment_cache.py    # Alignment cache
│       ├── aggregate_results.py          # Results aggregation
│       └── detect_primer_offsets.py      # UMI detection
│
└── archive/                               # Obsolete versions (for reference)
    ├── obsolete_workflow_scripts/        # 5 old script versions
    └── obsolete_snakefiles/              # 1 old Snakefile
```

---

## Verification

### ✅ No Version Suffixes Remain

```bash
$ find . -name "*_v[0-9]*" -type f | grep -v archive | grep -v __pycache__
# (no output - all cleaned up!)
```

### ✅ Current Workflow Scripts

```bash
$ ls -1 workflow/scripts/*.py | grep -v __pycache__
align_unique_sequences.py
aggregate_results.py
classify_per_guide_donor.py
collapse_reads.py                  # ← Renamed from v3
detect_primer_offsets.py
expand_to_samples.py               # ← Renamed from v2
extract_alignment_cache.py
extract_global_unique_seqs.py
```

### ✅ Current Snakefiles

```bash
$ ls -1 workflow/Snakefile*
Snakefile
Snakefile_cached                   # ← Renamed from v2
```

---

## Workflow Components

### Two-Level Caching System

**Level 1 (Alignment):** Each unique sequence aligned once
- Shared across ALL samples (reference constant)
- Cross-primer-pair sharing (KR2476 + KR2478)
- **Script:** `align_unique_sequences.py`

**Level 2 (Classification):** Each (sequence, guide, donor) classified once
- Reuses Level 1 alignment cache
- Only ~5 unique guide sequences, ~15 unique donor templates
- **Script:** `classify_per_guide_donor.py`

**Expansion:** Classifications mapped back to per-sample results
- **Script:** `expand_to_samples.py` (renamed from `expand_to_samples_v2.py`)

---

## Benefits

### 1. Clean File Naming
- No more `_v2`, `_v3` suffixes
- Single correct version of all files
- Clear workflow structure

### 2. Consistent References
- All projects updated to use clean names
- No confusion about which version to use
- Easy to maintain

### 3. Preserved History
- Obsolete versions archived for reference
- Can compare if needed
- Clean separation of current vs historical

---

## Testing Recommendations

### Test Workflow Execution

```bash
cd /oak/stanford/groups/larsms/Users/kevinroy/projects/HEK293/2025_HEK293_LGC_reporter

# Test workflow with renamed files
sbatch scripts/submit_trace_workflow.sh

# Should use:
# - Snakefile_cached (renamed from v2)
# - collapse_reads.py (renamed from v3)
# - expand_to_samples.py (renamed from v2)
```

### Verify Script Execution

```bash
# Test individual scripts still work
python3 /oak/.../trace/workflow/scripts/collapse_reads.py --help
python3 /oak/.../trace/workflow/scripts/expand_to_samples.py --help
```

---

## Archive Contents

### Obsolete Workflow Scripts (5 files, 48 KB)

```
archive/obsolete_workflow_scripts/
├── classify_sample_v2.py (10 KB)      - Unused versioned script
├── collapse_reads.py (8.2 KB)         - Old version
├── collapse_reads_v2.py (7.0 KB)      - Intermediate version
├── expand_to_samples.py (4.2 KB)      - Old version
└── expand_to_samples_v2.py.backup (4.7 KB) - Backup file
```

### Obsolete Snakefiles (1 file, 7.4 KB)

```
archive/obsolete_snakefiles/
└── Snakefile_cached (7.4 KB)          - Old cached workflow
```

---

## Related Documentation

- **Core Bug Fixes:** `/oak/.../projects/HEK293/.../TRACE_BUG_FIXES_SUMMARY.md`
- **Project Cleanup:** `/oak/.../projects/HEK293/.../CLEANUP_COMPLETE.md`
- **Version Cleanup:** `/oak/.../projects/HEK293/.../VERSION_CLEANUP_COMPLETE.md`
- **Workflow Guide:** `/oak/.../projects/HEK293/.../TRACE_WORKFLOW.md`

---

## Migration Notes

### If You Need to Revert

Old versions are preserved in `archive/` directories:

```bash
# Restore old version if needed
cp archive/obsolete_workflow_scripts/collapse_reads_v2.py workflow/scripts/
```

### If You're Using TRACE in Other Projects

Update your scripts to reference clean names:

```bash
# Find and replace in your project
sed -i 's|Snakefile_cached_v2|Snakefile_cached|g' *.sh
sed -i 's|collapse_reads_v3|collapse_reads|g' *.sh
sed -i 's|expand_to_samples_v2|expand_to_samples|g' *.sh
```

---

## Verification Checklist

- ✅ All `_v2`/`_v3` suffixes removed from active files
- ✅ Obsolete versions archived (not deleted)
- ✅ Snakefile_cached references updated scripts
- ✅ HEK293 project scripts updated (8 files)
- ✅ No remaining version suffixes in active code
- ✅ Archive directories created and populated
- ✅ File history preserved for reference

---

**Cleanup completed successfully. TRACE package now has clean, version-suffix-free naming.**
