# Obsolete Code Archive

This directory contains code that has been superseded by more efficient implementations.

## Files

### `kmer.py` - K-mer based HDR detection (OBSOLETE)
**Replaced by:** `core/edit_distance_hdr.py`

**Why obsolete:**
- K-mer approach required pre-generating all expected k-mers from donor template
- Performance issues: O(n²) string scanning for each k-mer
- Limited to SNV-only donors (no indel support)
- Superseded by edit-distance approach which:
  - Handles any number of SNVs without combinatorial explosion
  - Now supports donor indels via Needleman-Wunsch alignment
  - More efficient: O(n) single-pass classification per read

### `combinatorial_hdr.py` - Combinatorial HDR variant generator (OBSOLETE)
**Replaced by:** Edit-distance approach (no enumeration needed)

**Why obsolete:**
- Generated all possible combinations of donor SNVs for classification
- Memory exhaustion risk with >10 SNVs (2^N combinations)
- Not used in current workflow (cached v2 uses edit_distance_hdr.py)
- Edit-distance approach doesn't need to enumerate variants

## Migration Notes

The current TRACE workflow (`Snakefile_cached_v2`) exclusively uses `edit_distance_hdr.py`. These files can be safely removed if no legacy pipelines depend on them.

**Date archived:** 2026-03-09
