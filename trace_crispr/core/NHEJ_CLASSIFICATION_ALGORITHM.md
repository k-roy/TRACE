# Cut-Site-Aware NHEJ Classification Algorithm

**Date**: 2026-03-08
**Implementation**: `trace_crispr/core/edit_distance_hdr.py`
**Key Functions**: `_get_indel_position_range()`, `_could_involve_cut_site()`

---

## Overview

This document describes the cut-site-aware NHEJ classification algorithm that distinguishes true NHEJ repair products from sequencing artifacts and other non-NHEJ indels.

**Key Principle**: NHEJ indels MUST involve the Cas9 cut site. An indel can only be classified as NHEJ if its microhomology range includes the cut site position.

---

## The Problem

### Previous Approach (WRONG)

```python
# Old approach: Accept any indel within ±10bp window
if abs(indel_pos - cut_site) <= 10:
    outcome = "NHEJ"  # WRONG!
```

**Issues**:
1. **Distance-based, not mechanistic**: Being "near" the cut site doesn't mean the indel resulted from NHEJ
2. **Ignores microhomology**: Doesn't account for tandem repeats that create alignment ambiguity
3. **False positives**: Classifies homopolymer sequencing errors as NHEJ

**Example of False Positive**:
```
Cut site: position 222
Reference: ...TTTTTTT... at position 268 (7T homopolymer)
Indel: 1bp deletion at position 268
Distance: 268 - 222 = 46bp

Old algorithm: 46 > 10 → NOT NHEJ (correct by accident)
But with 10bp window: Would still be wrong for closer homopolymers
```

### Correct Approach (NEW)

```python
# New approach: Check if indel COULD involve cut site
min_pos, max_pos = get_microhomology_range(indel)
if (min_pos <= cut_site + 1) and (max_pos >= cut_site - 1):
    outcome = "NHEJ_INDEL"  # Must involve cut site!
```

**Key Insight**: The same biological indel can be represented at different positions in the alignment due to microhomology. We need to check if ANY of these possible positions includes the cut site.

---

## Algorithm Details

### Function 1: `_get_indel_position_range()`

**Purpose**: Calculate the range of positions where an indel could be placed due to microhomology-mediated alignment ambiguity.

#### For Deletions

```python
def _get_indel_position_range_deletion(indel_pos, indel_length, ref_seq):
    """
    Example:
    Reference: ATCGATCGATCG
    Deletion: ATC at position 3 (deletes "GAT")

    Can this deletion be shifted?
    - Left: Check if "GAT" prefix matches sequence before deletion
    - Right: Check if "GAT" suffix matches sequence after deletion

    If there's microhomology, the deletion could be at multiple positions.
    """
    deleted_seq = ref_seq[indel_pos:indel_pos + indel_length]

    # Check left shift potential
    left_shift = 0
    for i in range(1, min(indel_length, indel_pos) + 1):
        if deleted_seq[:i] == ref_seq[indel_pos - i:indel_pos]:
            left_shift = i
        else:
            break

    # Check right shift potential
    right_shift = 0
    for i in range(1, min(indel_length, len(ref_seq) - indel_pos - indel_length) + 1):
        if deleted_seq[-i:] == ref_seq[indel_pos + indel_length:indel_pos + indel_length + i]:
            right_shift = i
        else:
            break

    # Return range
    return (indel_pos - left_shift, indel_pos + indel_length + right_shift)
```

**Example**:
```
Reference:  ATCGATCGATCGATCG
            012345678901234567

Deletion at pos 3, length 3 (deletes "GAT"):
- Deleted sequence: "GAT"
- Can shift left? "GAT"[:1] = "G" vs ref[2:3] = "G" ✓ → left_shift = 1
- Can shift left more? "GAT"[:2] = "GA" vs ref[1:3] = "CG" ✗ → stop
- Can shift right? "GAT"[-1:] = "T" vs ref[6:7] = "T" ✓ → right_shift = 1
- Can shift right more? "GAT"[-2:] = "AT" vs ref[6:8] = "TC" ✗ → stop

Result: Deletion could be at positions [2, 7]
- At pos 2: delete "CGA"
- At pos 3: delete "GAT" (original alignment)
- At pos 4: delete "ATC"
- ...
```

#### For Insertions

```python
def _get_indel_position_range_insertion(indel_pos, indel_seq, ref_seq):
    """
    Accounts for tandem repeats:
    - Homopolymers: AAAAA
    - Dinucleotide repeats: ATATAT
    - Trinucleotide repeats: ATCATCATC
    - etc.

    Checks up to 5x insertion length for repeat context.
    """
    ins_len = len(indel_seq)

    # Check left extension
    left_extend = 0
    for shift in range(1, min(ins_len * 5, indel_pos) + 1):
        ref_context = ref_seq[indel_pos - shift:indel_pos]
        if ref_context == indel_seq[-shift:]:
            left_extend = shift
        else:
            break

    # Check right extension
    right_extend = 0
    for shift in range(1, min(ins_len * 5, len(ref_seq) - indel_pos) + 1):
        ref_context = ref_seq[indel_pos:indel_pos + shift]
        if ref_context == indel_seq[:shift]:
            right_extend = shift
        else:
            break

    return (indel_pos - left_extend, indel_pos + right_extend)
```

**Example (Homopolymer)**:
```
Reference: ...AAAAAAAAAA... (10 A's at position 100-109)
Insertion: "A" at position 105

Can slide within homopolymer:
- Left: A's extend from 100-104 (5 positions left)
- Right: A's extend from 105-109 (5 positions right)

Result: Insertion could be at [100, 110]
```

**Example (Dinucleotide Repeat)**:
```
Reference: ...ATATATATATAT... (6x AT repeat at position 200-211)
Insertion: "AT" at position 206

Can slide within repeat:
- Left: AT pattern extends 3 repeats left (6bp)
- Right: AT pattern extends 3 repeats right (6bp)

Result: Insertion could be at [200, 212]
```

### Function 2: `_could_involve_cut_site()`

**Purpose**: Check if an indel could plausibly result from NHEJ repair at the cut site.

```python
def _could_involve_cut_site(indel_pos, indel_length, indel_type,
                            indel_seq, ref_seq, cut_site, window=1):
    """
    NHEJ indels MUST involve the cut site.

    Args:
        window: Tiny window for Cas9 cutting position uncertainty (±1bp).
                SpCas9 cuts 3bp upstream of PAM, but can vary ±1bp.

    Returns:
        True if indel's microhomology range includes cut_site ± window
    """
    # Get the range where this indel could be placed
    min_pos, max_pos = _get_indel_position_range(
        indel_pos, indel_length, indel_type, indel_seq, ref_seq
    )

    # Check if cut site falls within range (with small window)
    return (min_pos <= cut_site + window) and (max_pos >= cut_site - window)
```

**Examples**:

```
Cut site: position 222

Example 1: Homopolymer indel (NOT NHEJ)
- Indel: 1bp deletion at position 268 (in 7T homopolymer)
- Range: [262, 274] (can slide within 7T run)
- Cut site 222 in [262, 274]? NO
- Result: NOT NHEJ (sequencing error)

Example 2: Near cut site but not involving (NOT NHEJ)
- Indel: 1bp deletion at position 220 (in 2T homopolymer)
- Range: [220, 221] (small tandem repeat)
- Cut site 222 in [220, 221]? NO
- Result: NOT NHEJ (too far, likely artifact)

Example 3: Involves cut site (IS NHEJ!)
- Indel: 2bp deletion at position 221
- Range: [221, 223] (deletion spans positions 221-222)
- Cut site 222 in [221, 223]? YES!
- Result: IS NHEJ (deletion involves the cut site)

Example 4: Large deletion spanning cut site (IS NHEJ!)
- Indel: 10bp deletion at position 218
- Range: [218, 228] (deletion spans positions 218-227)
- Cut site 222 in [218, 228]? YES!
- Result: IS NHEJ (deletion includes cut site)
```

---

## Integration with Classification

### Updated Outcome Categories

```python
def classify_read_edit_distance(..., nhej_quantification_window=1):
    # ... parse edits ...

    # Separate indels by whether they could involve cut site
    nhej_indels = []
    non_nhej_indels = []

    for edit in edits:
        if edit.edit_type != EditType.SNV and not edit.is_donor_encoded:
            indel_seq = edit.read_base if edit.edit_type == EditType.INSERTION else edit.ref_base

            could_be_nhej = _could_involve_cut_site(
                indel_pos=edit.ref_position,
                indel_length=edit.size,
                indel_type=edit.edit_type,
                indel_seq=indel_seq,
                ref_seq=ref_seq,
                cut_site=cut_site,
                window=nhej_quantification_window
            )

            if could_be_nhej:
                nhej_indels.append(edit)
            else:
                non_nhej_indels.append(edit)

    # Determine outcome
    if n_donor > 0 and len(nhej_indels) > 0:
        outcome = 'HDR_PLUS_NHEJ_INDEL'
    elif len(nhej_indels) > 0 and n_donor == 0:
        outcome = 'NHEJ_INDEL'
    elif len(non_nhej_indels) > 0 and n_donor == 0:
        outcome = 'NON_DONOR_NON_NHEJ_INDEL'
    # ... etc ...
```

### New Outcome Categories

| Category | Description | Example |
|----------|-------------|---------|
| **NHEJ_INDEL** | Indels that INVOLVE cut site | 2bp deletion spanning position 222 |
| **HDR_PLUS_NHEJ_INDEL** | Donor SNVs + NHEJ indels | HDR with concurrent NHEJ at cut site |
| **NON_DONOR_NON_NHEJ_INDEL** | Indels NOT at cut site | Homopolymer indel at +46bp from cut |
| **HDR_PLUS_OTHER** | Donor SNVs + non-NHEJ edits | HDR with sequencing errors |

---

## Parameters

### `nhej_quantification_window` (default: 1bp)

**Purpose**: Account for Cas9 cutting position uncertainty.

- SpCas9 cuts 3bp upstream of PAM
- Cut position can vary ±1bp depending on exact cleavage site
- Window of ±1bp allows for this biological uncertainty

**Example**:
```
PAM: NGG at position 225
Expected cut: position 222 (3bp upstream)
Actual cut: could be 221, 222, or 223

With window=1:
- Accept indels with range overlapping [221, 223]
```

### `snv_distance_filter` (default: 50bp)

**Purpose**: Only count non-donor SNVs within ±50bp of cut site.

**Rationale**: SNVs >50bp from cut site are likely sequencing errors, not editing outcomes.

### `homopolymer_filter` (default: 0, disabled)

**Purpose**: Filter 1bp indels at long homopolymer runs.

**Recommended**: 5-7 for Illumina data (filter 1bp indels at ≥5bp homopolymer runs).

**Example**:
```python
classify_read_edit_distance(
    ...,
    homopolymer_filter=6,  # Ignore 1bp indels at ≥6bp homopolymer runs
    nhej_quantification_window=1
)
```

---

## Testing

### Unit Test Examples

```python
def test_could_involve_cut_site():
    ref = "ATCGATCGATCGATCG"
    cut_site = 8  # Middle of sequence

    # Test 1: Deletion spanning cut site → NHEJ
    assert _could_involve_cut_site(
        indel_pos=7, indel_length=3, indel_type=EditType.DELETION,
        indel_seq="TCG", ref_seq=ref, cut_site=cut_site, window=1
    ) == True  # Range [7, 10] includes cut site 8

    # Test 2: Deletion far from cut site → NOT NHEJ
    assert _could_involve_cut_site(
        indel_pos=1, indel_length=2, indel_type=EditType.DELETION,
        indel_seq="TC", ref_seq=ref, cut_site=cut_site, window=1
    ) == False  # Range [1, 3] does NOT include cut site 8

    # Test 3: Insertion in tandem repeat → depends on range
    ref_with_repeat = "ATATATATATATAT"
    assert _could_involve_cut_site(
        indel_pos=6, indel_length=2, indel_type=EditType.INSERTION,
        indel_seq="AT", ref_seq=ref_with_repeat, cut_site=6, window=1
    ) == True  # Can slide within repeat to include cut site
```

---

## Performance Impact

### Computational Cost

- **Negligible**: Microhomology analysis adds ~1-5ms per indel
- **Benefit**: More accurate classification, fewer false positives

### Memory

- **No increase**: Algorithm operates on individual edits, no additional storage

---

## Biological Validation

### Control Samples

Expected outcome in negative controls (no Cas9 or no guide):
- **NHEJ_INDEL**: ~0% (no DSBs → no NHEJ)
- **NON_DONOR_NON_NHEJ_INDEL**: ~0.5-1% (sequencing errors, homopolymer artifacts)

### Edited Samples

Expected outcome in Cas9-edited samples:
- **NHEJ_INDEL**: 5-50% depending on repair efficiency
- **NON_DONOR_NON_NHEJ_INDEL**: ~0.5-1% (background errors)

**Key Insight**: If control samples show high NHEJ_INDEL rates, the algorithm is too permissive. With cut-site-aware classification, controls should show near-zero NHEJ_INDEL.

---

## References

### Microhomology-Mediated End Joining (MMEJ)

- Microhomology at deletion boundaries can cause alignment ambiguity
- Same biological deletion can be represented at multiple positions
- Affects indel calling in tandem repeat regions

### Cas9 Cutting Mechanism

- SpCas9 cuts 3bp upstream of PAM on target strand
- Blunt-ended DSB (some 5' overhang possible)
- NHEJ repair can result in indels at or near the cut site

### Illumina Sequencing Errors

- Homopolymer regions prone to 1bp insertion/deletion errors
- Error rate increases with homopolymer length
- Must distinguish from true editing events

---

## Future Improvements

### 1. Variable Window by Edit Size

Currently uses fixed ±1bp window. Could adapt based on indel size:
```python
window = min(1 + indel_length // 10, 3)  # Larger indels → larger window
```

### 2. Sequence Context Scoring

Weight indels by sequence context complexity:
- High complexity (unique sequence) → higher confidence
- Low complexity (tandem repeat) → lower confidence

### 3. Machine Learning Classification

Train classifier on validated NHEJ vs artifact indels using features:
- Distance to cut site
- Microhomology range
- Sequence complexity
- Indel size
- Base quality scores

---

## Changelog

### 2026-03-08: Initial Implementation
- Added `_get_indel_position_range()` for microhomology analysis
- Added `_could_involve_cut_site()` for cut-site-aware filtering
- Updated outcome categories to NHEJ_INDEL, NON_DONOR_NON_NHEJ_INDEL, etc.
- Changed default `nhej_quantification_window` from 10bp to 1bp
- Added tandem repeat detection for insertions

---

For implementation details, see:
- `trace_crispr/core/edit_distance_hdr.py:186-333`
- `workflow/scripts/classify_per_guide_donor.py` (integration)
