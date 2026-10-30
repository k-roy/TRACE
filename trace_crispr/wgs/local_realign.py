"""Local realignment of the reads at a localized edit site to its wild-type and design haplotypes.

Why this exists (SGTC planning/328, 2026-10-30). bwa soft-clips a read across a dense design: ten substitutions
in 30 bp score worse than a clip. An allele reader that reads only aligned bases then sees NO evidence at the
designed bases of a correctly edited clone. On the NRD1/NAB3/SEN1 NNS cohort, 119 of 139 correct edits the
aligned-base reader could not call carried the design in clipped bases, and 1,090 of the 1,198 clips at the
donor windows sat on the designed span.

What it does. Once the edit site is known (its intended alleles), every primary read whose aligned span OR
soft clips reach the designed span is realigned, with its FULL query sequence, to two local haplotypes: the
reference, and the reference carrying every intended allele. Scoring is bwa-mem's defaults (match 1, mismatch 4,
gap open 6, gap extend 1, clipping 5 per end), so any change from the genome-wide alignment comes from the
haplotype, not from a different scoring scheme, and a read can still clip: a plasmid-backbone or adapter junction
stays clipped instead of being forced into an allele call.

How an allele is read. A substitution allele is read by IDENTITY (alt, ref or other) from the better-scoring of
the two alignments, never by which haplotype won, so partial incorporation (some designed bases installed, others
wild type) stays visible. A length-changing allele is uninformative for a read whose score does not change when
that allele alone is removed (inside a repeat it does not span); otherwise it too is read by exact identity across
the allele in the better alignment.

Two guards keep the realignment from inventing evidence. (1) ANCHOR: a realignment is used only if every base the
genome-wide aligner placed outside the designed span (more than `ANCHOR_MARGIN` bp from it) keeps its reference
position, and at least `MIN_ANCHOR` such bases do; otherwise the read was moved, typically by sliding a nearby
error onto a designed base inside a repeat, and its realignment is discarded. (2) FALLBACK: where the realignment
is discarded, or reads an allele as uncovered or uninformative, the caller's aligned-base observation for that
read and allele is used instead, so every read the aligned-base reader could read is still read. The realignment
therefore only ever ADDS evidence the alignment hid, or replaces a reading with a better-scoring anchored one.

All coordinates are 1-based inclusive at the interface, as in `reconcile`.
"""
from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple

MATCH, MISMATCH, GAP_OPEN, GAP_EXTEND, CLIP = 1, 4, 6, 1, 5
#: Reference bases either side of a read's projected span searched by the realignment.
BAND = 40
#: Reads are fetched this far either side of the designed span, so a read aligned entirely outside it whose
#: soft clip reaches into it is still found.
FETCH_PAD = 300
#: Anchor rule: originally aligned bases farther than this from the designed span must keep their position.
ANCHOR_MARGIN = 5
MIN_ANCHOR = 5
_NEG = -10 ** 9


def align(query: str, target: str, expect_end: Optional[int] = None) -> Tuple[int, List[Tuple[int, Optional[int]]]]:
    """Clip-allowed local alignment of `query` to `target` with free target ends (Gotoh, bwa-mem scores).

    Returns `(score, columns)`, where `columns` lists `(query_index, target_index or None)` for every aligned
    query base in order; `None` marks an inserted query base. Clipped query bases do not appear; deleted
    target bases are implied by gaps in the target indices.

    `expect_end` is the target index where the read's original placement puts its last base. Among end cells
    of equal score it picks the one nearest that placement, so in a homopolymer or tandem repeat, where many
    placements score the same, realignment keeps the genome-wide aligner's placement instead of an arbitrary
    one.
    """
    n, m = len(query), len(target)
    H = [[_NEG] * (m + 1) for _ in range(n + 1)]
    E = [[_NEG] * (m + 1) for _ in range(n + 1)]  # gap in query (target base deleted), horizontal move
    F = [[_NEG] * (m + 1) for _ in range(n + 1)]  # gap in target (query base inserted), vertical move
    Hp = [[3] * (m + 1) for _ in range(n + 1)]     # 0 diag, 1 from E, 2 from F, 3 begin (after a clip)
    Ep = [[0] * (m + 1) for _ in range(n + 1)]     # 0 open from H, 1 extend
    Fp = [[0] * (m + 1) for _ in range(n + 1)]
    go_ge = GAP_OPEN + GAP_EXTEND
    for i in range(n + 1):
        begin = 0 if i == 0 else -CLIP
        Hi, Ei, Fi, Hpi, Epi, Fpi = H[i], E[i], F[i], Hp[i], Ep[i], Fp[i]
        Him1 = H[i - 1] if i else None
        Fim1 = F[i - 1] if i else None
        qc = query[i - 1] if i else ""
        for j in range(m + 1):
            best, ptr = begin, 3
            if j:
                e_open, e_ext = Hi[j - 1] - go_ge, Ei[j - 1] - GAP_EXTEND
                if e_ext > e_open:
                    Ei[j], Epi[j] = e_ext, 1
                else:
                    Ei[j], Epi[j] = e_open, 0
                if Ei[j] > best:
                    best, ptr = Ei[j], 1
            if i:
                f_open, f_ext = Him1[j] - go_ge, Fim1[j] - GAP_EXTEND
                if f_ext > f_open:
                    Fi[j], Fpi[j] = f_ext, 1
                else:
                    Fi[j], Fpi[j] = f_open, 0
                if Fi[j] > best:
                    best, ptr = Fi[j], 2
                if j:
                    d = Him1[j - 1] + (MATCH if qc == target[j - 1] and qc != "N" else -MISMATCH)
                    if d > best:
                        best, ptr = d, 0
            Hi[j], Hpi[j] = best, ptr
    best, bi, bj, bdist = _NEG, 0, 0, 0
    for i in range(n + 1):
        tail = 0 if i == n else -CLIP
        Hi = H[i]
        for j in range(m + 1):
            v = Hi[j] + tail
            if v < best:
                continue
            dist = abs(j + (n - i) - 1 - expect_end) if expect_end is not None else 0
            if v > best or dist < bdist:
                best, bi, bj, bdist = v, i, j, dist
    cols: List[Tuple[int, Optional[int]]] = []
    i, j, state = bi, bj, "H"
    while i > 0 or j > 0:
        if state == "H":
            p = Hp[i][j]
            if p == 3:
                break
            if p == 0:
                cols.append((i - 1, j - 1))
                i, j = i - 1, j - 1
            elif p == 1:
                state = "E"
            else:
                state = "F"
        elif state == "E":
            p = Ep[i][j]
            j -= 1
            if p == 0:
                state = "H"
        else:
            p = Fp[i][j]
            cols.append((i - 1, None))
            i -= 1
            if p == 0:
                state = "H"
    cols.reverse()
    return best, cols


class _Haplotypes:
    """The reference window and the window carrying every intended allele, with allele spans in each."""

    def __init__(self, fasta, chrom: str, alleles: Sequence[Tuple[int, str, str]]):
        self.alleles = sorted(alleles)
        lo = min(p for p, _, _ in self.alleles)
        hi = max(p + len(r) - 1 for p, r, _ in self.alleles)
        clen = fasta.get_reference_length(chrom)
        self.start = max(1, lo - FETCH_PAD - BAND)               # 1-based position of index 0
        self.end = min(clen, hi + FETCH_PAD + BAND)
        self.span = (lo, hi)
        self.wt = fasta.fetch(chrom, self.start - 1, self.end).upper()
        parts, cursor, shift = [], self.start, 0
        self.wt_spans, self.design_spans = [], []
        for pos, ref, alt in self.alleles:
            observed = self.wt[pos - self.start: pos - self.start + len(ref)]
            if observed != ref:
                raise ValueError(f"{chrom}:{pos}: intended REF {ref!r} does not match reference {observed!r}")
            parts.append(self.wt[cursor - self.start: pos - self.start])
            a = pos - self.start
            self.wt_spans.append((a, a + len(ref)))
            self.design_spans.append((a + shift, a + shift + len(alt)))
            parts.append(alt)
            shift += len(alt) - len(ref)
            cursor = pos + len(ref)
        parts.append(self.wt[cursor - self.start:])
        self.design = "".join(parts)
        self._shifts = [(p, len(a) - len(r)) for p, r, a in self.alleles]
        # reference position (1-based) of every design-haplotype base; None for an inserted base
        self.design_ref: List[Optional[int]] = []
        cursor = self.start
        for pos, ref, alt in self.alleles:
            self.design_ref.extend(range(cursor, pos))
            if len(ref) == len(alt):
                self.design_ref.extend(range(pos, pos + len(ref)))
            else:
                self.design_ref.extend([pos] + [None] * (len(alt) - 1))
            cursor = pos + len(ref)
        self.design_ref.extend(range(cursor, self.end + 1))
        assert len(self.design_ref) == len(self.design)

    def design_index(self, pos1: int) -> int:
        """Design-haplotype index of reference position `pos1` (for an unchanged base)."""
        return pos1 - self.start + sum(d for p, d in self._shifts if p < pos1)

    def single_allele(self, k: int, base: str) -> str:
        """`base` haplotype ("wt" or "design") with allele k flipped to the other state."""
        pos, ref, alt = self.alleles[k]
        if base == "wt":
            a, b = self.wt_spans[k]
            return self.wt[:a] + alt + self.wt[b:]
        a, b = self.design_spans[k]
        return self.design[:a] + ref + self.design[b:]


def _read_span(read) -> Tuple[int, int]:
    """1-based projected span of a read: its aligned span extended by its soft clips."""
    ops = read.cigartuples or []
    lead = ops[0][1] if ops and ops[0][0] == 4 else 0
    trail = ops[-1][1] if ops and ops[-1][0] == 4 else 0
    return read.reference_start + 1 - lead, read.reference_end + trail


def _observe(cols, span, seq, quals, min_base_qual) -> Optional[str]:
    """Read bases aligned to target `span` [a, b), with insertions after any base in it; None if not covered."""
    a, b = span
    tpos = [t for _, t in cols if t is not None]
    if not tpos or min(tpos) > a or max(tpos) < b - 1:
        return None
    out, inside, last_t = [], False, None
    for q, t in cols:
        if t is not None:
            inside = a <= t < b
            last_t = t
            if inside:
                if quals[q] < min_base_qual or seq[q] not in "ACGT":
                    return None
                out.append(seq[q])
        elif last_t is not None and a <= last_t < b:
            if quals[q] < min_base_qual or seq[q] not in "ACGT":
                return None
            out.append(seq[q])
    return "".join(out)


def _anchored(read, cols, to_ref, lo: int, hi: int) -> bool:
    """True if the realignment kept every originally aligned base outside the designed span in place."""
    new = {q: to_ref(t) for q, t in cols if t is not None}
    kept = 0
    for q, r0 in read.get_aligned_pairs(matches_only=True):
        r1 = r0 + 1
        if lo - ANCHOR_MARGIN <= r1 <= hi + ANCHOR_MARGIN or q not in new:
            continue
        if new[q] != r1:
            return False
        kept += 1
    return kept >= MIN_ANCHOR


def realign_read(read, haps: "_Haplotypes", min_base_qual: int) -> Tuple[bool, List[Optional[str]]]:
    """Realign one read to both haplotypes; `(anchored, observation per allele in haps.alleles order)`."""
    s, e = _read_span(read)
    seq, quals = read.query_sequence.upper(), read.query_qualities
    w0 = max(0, s - BAND - haps.start)
    w1 = min(len(haps.wt), e + BAND - haps.start + 1)
    d0 = max(0, haps.design_index(max(haps.start, s - BAND)))
    d1 = min(len(haps.design), haps.design_index(min(haps.end, e + BAND)) + 1)
    wt_end, ds_end = e - haps.start - w0, haps.design_index(min(e, haps.end)) - d0
    wt_score, wt_cols = align(seq, haps.wt[w0:w1], expect_end=wt_end)
    ds_score, ds_cols = align(seq, haps.design[d0:d1], expect_end=ds_end)
    use_design = ds_score > wt_score
    if use_design:
        cols = [(q, None if t is None else t + d0) for q, t in ds_cols]
        to_ref = lambda t: haps.design_ref[t]
    else:
        cols = [(q, None if t is None else t + w0) for q, t in wt_cols]
        to_ref = lambda t: t + haps.start
    lo, hi = haps.span
    obs: List[Optional[str]] = [None] * len(haps.alleles)
    if not _anchored(read, cols, to_ref, lo, hi):
        return False, obs
    here = ds_score if use_design else wt_score
    for k, (pos, ref, alt) in enumerate(haps.alleles):
        span = (haps.design_spans if use_design else haps.wt_spans)[k]
        if len(ref) == len(alt):
            o = _observe(cols, span, seq, quals, min_base_qual)
            if o is not None:
                obs[k] = "alt" if o == alt else "ref" if o == ref else "other"
            continue
        if _observe(cols, (span[0], span[1] + 1), seq, quals, min_base_qual) is None:
            continue
        flipped = haps.single_allele(k, "design" if use_design else "wt")
        delta = (len(ref) - len(alt)) if use_design else (len(alt) - len(ref))
        lo_i, hi_i = (d0, d1) if use_design else (w0, w1)
        other_score, _ = align(seq, flipped[lo_i: hi_i + delta], expect_end=(ds_end if use_design else wt_end))
        if other_score == here:
            continue                      # the read cannot tell this allele from its absence (e.g. inside a repeat)
        # Exact identity across the allele in the chosen alignment, as the aligned-base reader requires: a read
        # whose anchor base differs from the declared one is `other`, not the indel. A repeat-shifted
        # representation of the same indel is the same haplotype and aligns to the design without a gap.
        o = _observe(cols, span, seq, quals, min_base_qual)
        if o is not None:
            obs[k] = "alt" if o == alt else "ref" if o == ref else "other"
    return True, obs


def realigned_allele_counts(bam, fasta, chrom: str, alleles: Sequence[Tuple[int, str, str]], *,
                            min_base_qual: int = 13, min_mapq: int = 20, fallback=None,
                            fragments_out: Optional[list] = None,
                            ) -> Tuple[List[Dict[str, int]], int, Dict[str, int]]:
    """Fragment counts per intended allele from local realignment, in the order of `alleles`.

    `fallback(read, k)` returns the aligned-base observation ("alt", "ref", "other" or None) of allele k for a
    read; it is used where the realignment is discarded (not anchored) or reads that allele as uncovered.
    Returns `(counts, n_unmapped_at_locus, read_tally)`. Each count has alt, ref, other, depth,
    low_mapq_reads, unknown_mapq_reads. Read filters match `reconcile.intended_allele_evidence`: primary,
    mapped, non-duplicate, QC-passing, MAPQ != 255 and >= min_mapq. Placed-unmapped reads in the fetch window
    are COUNTED, not used: without an alignment of their own they are not local evidence. `read_tally` counts
    reads by how they were read (anchored realignment, displaced and so fallen back, observations filled from
    the fallback). If `fragments_out` is a list, one observation list per merged fragment (in `alleles` order:
    "alt", "ref", another base, "conflict" or None) is appended to it, for haplotype-level calling.
    """
    order = sorted(range(len(alleles)), key=lambda k: tuple(alleles[k]))
    haps = _Haplotypes(fasta, chrom, [tuple(alleles[k]) for k in order])
    lo, hi = haps.span
    n = len(alleles)
    fragments: Dict[tuple, List[Optional[str]]] = {}
    low = unknown = unmapped = 0
    tally = {"reads_realigned_anchored": 0, "reads_realign_displaced": 0, "observations_from_fallback": 0,
             "observations_from_realignment": 0}
    for index, read in enumerate(bam.fetch(chrom, max(0, lo - 1 - FETCH_PAD), hi + FETCH_PAD)):
        if read.flag & 0x4:
            unmapped += 1
            continue
        if read.flag & 0xF00 or read.query_sequence is None or read.query_qualities is None:
            continue
        s, e = _read_span(read)
        if e < lo or s > hi:
            continue
        if read.mapping_quality == 255:
            unknown += 1
            continue
        if read.mapping_quality < min_mapq:
            low += 1
            continue
        anchored, sorted_obs = realign_read(read, haps, min_base_qual)
        tally["reads_realigned_anchored" if anchored else "reads_realign_displaced"] += 1
        obs: List[Optional[str]] = [None] * n
        for k_sorted, k in enumerate(order):
            o = sorted_obs[k_sorted]
            if o is not None:
                tally["observations_from_realignment"] += 1
            elif fallback is not None:
                o = fallback(read, k)
                if o is not None:
                    tally["observations_from_fallback"] += 1
            obs[k] = o
        if all(o is None for o in obs):
            continue
        key = (read.get_tag("RG") if read.has_tag("RG") else "", read.query_name) if read.is_paired else (None, index)
        if key in fragments:
            prev = fragments[key]
            for k in range(n):
                if obs[k] is None:
                    continue
                if prev[k] is None:
                    prev[k] = obs[k]
                elif prev[k] != obs[k]:
                    prev[k] = "conflict"
        else:
            fragments[key] = obs
    if fragments_out is not None:
        fragments_out.extend(list(f) for f in fragments.values())
    counts = []
    for k in range(n):
        values = [f[k] for f in fragments.values() if f[k] not in (None, "conflict")]
        alt_n, ref_n = values.count("alt"), values.count("ref")
        counts.append({"alt": alt_n, "ref": ref_n, "other": len(values) - alt_n - ref_n, "depth": len(values),
                       "low_mapq_reads": low, "unknown_mapq_reads": unknown})
    return counts, unmapped, tally
