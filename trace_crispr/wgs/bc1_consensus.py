"""
trace_crispr.wgs.bc1_consensus
==============================

Call the **bc1 lineage barcode(s)** of a pooled-MAGESTIC clone directly from its
whole-genome BAM, so each cell can be bridged to its designed donor via the
``bc1 -> donor`` reference table.

Design (2026-06-04) — *not* a port of Shengdi Li's trie caller
-------------------------------------------------------------
The barcode is a 20 bp degenerate locus flanked by two long **constant** junction
sequences, on an integrated contig (and an episomal-plasmid copy). With 75 bp Tn5
reads a single read captures the barcode plus part of **one** flank, never both,
and coverage is near-singleton (most calls rest on a single read). So multi-read
consensus machinery (Shengdi's prefix trie + ``Tscore``/``Uscore``) is moot and
buggy at this depth. Verified on the real BAMs: the aligner **soft-clips** the
random barcode (CIGARs like ``3S38M34S`` — ~38 bp anchored in the flank, the
barcode in the soft-clip), so a reference-column pileup would recover nothing.

This module instead does robust, deterministic, **per-read** work:

1. anchor on the **constant flank** (never align against the random barcode);
2. map the flank↔barcode boundary to a query offset via the CIGAR and slice the
   20 bp directly by query index, **recovering soft-clipped barcode bases**;
3. drop read-2 amplicon-hop contaminants (mate off the barcode contig);
4. dedup, then bridge the extracted 20-mer to the designed ``bc1 -> donor`` set.

Confidence is feature-based (bridge tier + anchor cleanliness + support), never a
read count. The genome/plasmid tracks are called independently and reported
separately (their flanks differ, so two-track agreement is genuine, not circular).
The double-edit-vs-contamination call and the edit-VAF corroboration live in a
**separate reconciliation step** (they need the variant calls), per Kevin's call.

Provenance: the constant-flank-anchor insight and the contaminant fingerprint are
from Shengdi Li's ``call_barcode_from_bam.pl`` / ``filter_barcode_region.pl``; the
extraction, bridge, and confidence are a ground-up redesign (Kevin R. Roy, 2026).
"""

from __future__ import annotations

import argparse
import glob
import os
import re
import sys
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import pysam

from ..core.cigar import parse_cigar_to_operations

# The genuine barcode read is a proper pair (flags 83/99/147/163, == read.is_proper_pair
# for a primary mapped read), which extract_track() requires; see is_amplicon_read /
# is_read2_hop for the contaminant gates layered on top.

_RC_TABLE = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq: str) -> str:
    """Reverse-complement (preserving N)."""
    return seq.translate(_RC_TABLE)[::-1]


@dataclass
class TrackConfig:
    """One barcode track (genomic-integrated or episomal-plasmid)."""

    name: str  # e.g. "genome" / "plasmid"
    contig: str
    bc_start: int  # 1-based inclusive (first barcode base)
    bc_end: int  # 1-based inclusive (last barcode base)
    left_junction: str  # constant left flank (this contig's forward orientation)
    right_junction: str  # constant right flank (this contig's forward orientation)
    revcomp: bool = False  # True if forward-on-contig barcode is the rev-comp of the
    # canonical bc1. Canonical (== the bc1->donor reference table / SHA345 / yT182 truth)
    # is revcomp(genome-forward): Shengdi's raw genome call is rc(truth) and is rev-comped
    # downstream. The plasmid cassette is inverted vs the genome integrant
    # (plasmid left flank == revcomp(genome right flank), verified 150 bp exact), so
    # plasmid-forward IS already canonical. => genome track revcomp=True, plasmid=False.
    # (Plasmid==canonical-vs-truth is inferred from the flank inversion, not yet observed
    # against a confident plasmid call.)

    # Optional bc1-amplicon contaminant telltale (Shengdi's filter_barcode_region.pl). The
    # bc1 fitness amplicon cross-contaminates WGS, including as PROPER PAIRS mapped on-contig
    # (a contaminant class the read-2-hop / mate-off-contig rule cannot catch). Such reads
    # have an aligned outer boundary at a fixed amplicon-primer coordinate plus a randomer
    # soft-clip/indel there; genuine Tn5 reads have random boundaries. If a window is set, a
    # read whose aligned start is in amplicon_left_window with a defect in the first
    # amplicon_left_defect query bases, OR whose aligned end is in amplicon_right_window with
    # a defect in the last amplicon_right_defect bases, is dropped. None => filter disabled.
    amplicon_left_window: Optional[Tuple[int, int]] = None  # 1-based inclusive start window
    amplicon_right_window: Optional[Tuple[int, int]] = None  # 1-based inclusive end window
    amplicon_left_defect: int = 5
    amplicon_right_defect: int = 7

    @property
    def bc_len(self) -> int:
        return self.bc_end - self.bc_start + 1

    def to_canonical(self, query_barcode: str) -> str:
        """Map a forward-on-contig barcode to the canonical bc1 orientation (the
        bc1->donor reference table / SHA345 / yT182 truth orientation). For a track
        with ``revcomp=True`` (the genome integrant) that is the reverse-complement of
        the forward-on-contig read; for ``revcomp=False`` (the plasmid) it is unchanged."""
        return revcomp(query_barcode) if self.revcomp else query_barcode

    def validate_against_reference(self, fasta: "pysam.FastaFile") -> None:
        """Assert the locus is an N-run flanked by the configured junctions.

        Fails LOUDLY on reference drift (the bug class that silently broke
        Shengdi's hard-coded coordinates).
        """
        if self.contig not in fasta.references:
            raise ValueError(f"track {self.name}: contig {self.contig!r} not in reference")
        # 0-based half-open fetch
        bc = fasta.fetch(self.contig, self.bc_start - 1, self.bc_end).upper()
        if set(bc) != {"N"}:
            raise ValueError(
                f"track {self.name}: {self.contig}:{self.bc_start}-{self.bc_end} is not an N-run ({bc!r})"
            )
        left = fasta.fetch(self.contig, self.bc_start - 1 - len(self.left_junction), self.bc_start - 1).upper()
        if left != self.left_junction.upper():
            raise ValueError(f"track {self.name}: left junction does not match reference")
        right = fasta.fetch(self.contig, self.bc_end, self.bc_end + len(self.right_junction)).upper()
        if right != self.right_junction.upper():
            raise ValueError(f"track {self.name}: right junction does not match reference")


# --------------------------------------------------------------------------- #
# Extraction kernel (the ~40 lines that matter)                                #
# --------------------------------------------------------------------------- #


def ref_to_query_offset(read: "pysam.AlignedSegment", ref0: int) -> Optional[int]:
    """Query index (0-based) aligned to — or **extrapolated to** — reference
    coordinate ``ref0`` (0-based), recovering soft-clipped barcode bases.

    If ``ref0`` falls inside a match op, return the exact aligned query index.
    If it falls in (or past) a soft-clip/insertion that abuts a match op, walk
    1:1 into the clip from the nearest match boundary — valid because the barcode
    is soft-clipped precisely *because* it diverges from the N-reference, so its
    bases sit contiguously in the query just past the alignment edge.

    Returns ``None`` if ``ref0`` cannot be related to any match op.
    """
    ops = parse_cigar_to_operations(read)
    matches = [op for op in ops if op.op_code in (0, 7, 8)]  # M, =, X
    if not matches:
        return None
    for op in matches:
        if op.ref_start <= ref0 < op.ref_end:  # inside an aligned block
            return op.query_start + (ref0 - op.ref_start)
    # extrapolate from the nearest match boundary in the direction of ref0
    first, last = matches[0], matches[-1]
    if ref0 >= last.ref_end:  # to the right of all matches -> trailing clip
        return last.query_end + (ref0 - last.ref_end)
    if ref0 < first.ref_start:  # to the left of all matches -> leading clip
        return first.query_start - (first.ref_start - ref0)
    return None  # falls in a deletion/intron between matches


@dataclass
class Extraction:
    """A barcode extracted from one read, with raw per-read features.

    ``barcode`` is the **canonical** bc1 (genome-forward; bridge/agreement use this).
    ``query_barcode`` is the raw forward-on-contig slice (== barcode for the genome
    track; == revcomp(barcode) for the inverted plasmid track) — kept for QC/debug.
    """

    barcode: str  # canonical bc1 orientation
    query_barcode: str  # raw forward-on-contig slice
    anchor: str  # "left" or "right"
    read_name: str
    is_reverse: bool
    mapq: int
    flank_mismatches: int  # mismatches in the anchored constant flank (frame QC)
    barcode_has_n: bool
    dedup_key: Tuple


def extract_barcode(read: "pysam.AlignedSegment", track: TrackConfig) -> Optional[Extraction]:
    """Extract the 20 bp barcode from one spanning read by flank-boundary query
    slicing. Returns ``None`` if the read does not anchor a flank cleanly.

    Routing: a read aligned starting left of the barcode anchors the LEFT flank
    (barcode in the trailing soft-clip); one aligned within/after the barcode
    anchors the RIGHT flank (barcode in the leading soft-clip).
    """
    if read.query_sequence is None or read.reference_start is None:
        return None
    bc_start0 = track.bc_start - 1  # 0-based first barcode base
    bc_end0 = track.bc_end  # 0-based one-past-last barcode base
    n = track.bc_len

    if read.reference_start < bc_start0:
        anchor = "left"
        q = ref_to_query_offset(read, bc_start0)  # first barcode base
        if q is None or q < 0:  # guard: a negative slice start wraps to the read 3' end
            return None
        barcode = read.query_sequence[q : q + n]
    else:
        anchor = "right"
        q_end = ref_to_query_offset(read, bc_end0)  # first base past barcode
        if q_end is None or q_end - n < 0:  # guard against a negative slice start
            return None
        barcode = read.query_sequence[q_end - n : q_end]

    if len(barcode) != n:
        return None  # read did not reach a full barcode width
    query_barcode = barcode.upper()
    canonical = track.to_canonical(query_barcode)

    flank_mm = _flank_mismatches(read, track, anchor)
    return Extraction(
        barcode=canonical,
        query_barcode=query_barcode,
        anchor=anchor,
        read_name=read.query_name,
        is_reverse=bool(read.is_reverse),
        mapq=read.mapping_quality,
        flank_mismatches=flank_mm,
        barcode_has_n=("N" in canonical),
        dedup_key=_dedup_key(read),
    )


def _flank_mismatches(read: "pysam.AlignedSegment", track: TrackConfig, anchor: str) -> int:
    """Count substitution mismatches in the aligned portion of the anchored
    constant flank — the main per-read frame-quality number."""
    if anchor == "left":
        region = (track.bc_start - 1 - len(track.left_junction), track.bc_start - 1)
        flank = track.left_junction.upper()
    else:
        region = (track.bc_end, track.bc_end + len(track.right_junction))
        flank = track.right_junction.upper()
    r0, r1 = region
    seq = read.query_sequence
    mm = 0
    for qpos, rpos in read.get_aligned_pairs(matches_only=True):
        if r0 <= rpos < r1:
            exp = flank[rpos - r0]
            got = seq[qpos].upper()
            if exp != "N" and got != "N" and exp != got:
                mm += 1
    return mm


def _dedup_key(read: "pysam.AlignedSegment") -> Tuple:
    """Collapse PCR duplicates to one *molecule*. Uses the fragment signature when
    the mate is mapped, else a single-read fallback (ref span + cigar). TLEN is
    abs()'d and read1/read2 are NOT distinguished, so both ends of the same fragment
    (e.g. a left- and a right-anchored mate that both extract the barcode) collapse
    to a single molecule rather than being counted twice."""
    if read.is_paired and not read.mate_is_unmapped and read.next_reference_id == read.reference_id:
        return ("pair", min(read.reference_start, read.next_reference_start), abs(read.template_length))
    return ("single", read.reference_start, read.reference_end, read.cigarstring)


# --------------------------------------------------------------------------- #
# Contaminant filter + per-sample extraction                                   #
# --------------------------------------------------------------------------- #


def is_read2_hop(read: "pysam.AlignedSegment", track: TrackConfig) -> bool:
    """Read-2 amplicon-hop contaminant test: a barcode-bearing read whose mate maps
    **off** the barcode contig. For a genuine integration read both ends sit near the
    locus (mate on-contig); a hopped bc1-amplicon read-2 has its genuine WGS read-1 mate
    on a real chromosome. NB this catches only the *cross-contig* hop; full bc1-amplicon
    fragments mapped as a proper pair ON the barcode contig evade it and require
    :func:`is_amplicon_read`."""
    if read.mate_is_unmapped:
        return True
    return read.next_reference_name != track.contig


def _boundary_defect(read: "pysam.AlignedSegment", n_bp: int, from_end: bool) -> bool:
    """True if a soft-clip/insertion/deletion/X op falls within the first (or, if
    ``from_end``, last) ``n_bp`` query bases — the amplicon-primer randomer signature
    (== Shengdi's 'defect in first 5bp / last 7bp' rule)."""
    ops = read.cigartuples or []
    if from_end:
        ops = list(reversed(ops))
    qpos = 0
    for op, length in ops:
        if op in (1, 2, 4, 8):  # I, D, S, X
            if qpos < n_bp:
                return True
        if op in (0, 1, 4, 7, 8):  # query-consuming: M, I, S, =, X
            qpos += length
        if qpos >= n_bp:
            break
    return False


def is_amplicon_read(read: "pysam.AlignedSegment", track: TrackConfig) -> bool:
    """bc1-amplicon contaminant telltale (Shengdi's ``filter_barcode_region.pl``):
    an aligned outer boundary at a fixed amplicon-primer coordinate **plus** a randomer
    soft-clip/indel at that boundary. Disabled (returns False) when no window is set."""
    if read.reference_start is None or read.reference_end is None:
        return False
    if track.amplicon_left_window is not None:
        lo, hi = track.amplicon_left_window
        start1 = read.reference_start + 1  # 1-based aligned start
        if lo <= start1 <= hi and _boundary_defect(read, track.amplicon_left_defect, from_end=False):
            return True
    if track.amplicon_right_window is not None:
        lo, hi = track.amplicon_right_window
        end1 = read.reference_end  # 1-based inclusive aligned end (pysam end is 0-based exclusive)
        if lo <= end1 <= hi and _boundary_defect(read, track.amplicon_right_defect, from_end=True):
            return True
    return False


@dataclass
class TrackCall:
    """All extractions for one (sample, track), pre-bridge."""

    track: str
    extractions: List[Extraction] = field(default_factory=list)
    n_reads_seen: int = 0
    n_contaminant: int = 0
    n_amplicon: int = 0  # subset of n_contaminant dropped by the amplicon telltale

    def barcode_support(self) -> Dict[str, int]:
        """Distinct-molecule (deduped) support per extracted barcode."""
        seen: Dict[str, set] = {}
        for e in self.extractions:
            seen.setdefault(e.barcode, set()).add(e.dedup_key)
        return {bc: len(keys) for bc, keys in seen.items()}


def extract_track(
    bam: "pysam.AlignmentFile",
    track: TrackConfig,
    window_pad: int = 160,
) -> TrackCall:
    """Stream the barcode window of one track and extract a barcode per spanning
    read, dropping read-2 hops. (pysam windowed fetch — efficient vs streaming
    the whole BAM.)"""
    call = TrackCall(track=track.name)
    bc0 = track.bc_start - 1  # 0-based first barcode base
    bc_end0 = track.bc_end  # 0-based one-past-last barcode base
    start = max(0, bc0 - window_pad)
    end = track.bc_end + window_pad
    for read in bam.fetch(track.contig, start, end):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if read.reference_start is None or read.reference_end is None:
            continue
        # A read can carry the barcode iff its SOFT-CLIP-EXTENDED reference span covers the
        # barcode. The barcode is soft-clipped (it diverges from the N-reference), so a
        # right-anchored read aligning the right flank starts AT bc_end with the barcode
        # entirely in its leading soft-clip — an aligned-span overlap test wrongly drops it.
        cig = read.cigartuples or []
        lead_s = cig[0][1] if cig and cig[0][0] == 4 else 0
        trail_s = cig[-1][1] if cig and cig[-1][0] == 4 else 0
        if read.reference_end + trail_s <= bc0 or read.reference_start - lead_s >= bc_end0:
            continue
        call.n_reads_seen += 1
        # Genuine barcode read = PROPER PAIR with mate on-contig. Drop (a) non-proper-pair
        # and cross-contig read-2 hops, and (b) on-contig bc1-amplicon proper pairs via the
        # amplicon-primer telltale (only active when the track has amplicon windows set).
        if not read.is_proper_pair or is_read2_hop(read, track):
            call.n_contaminant += 1
            continue
        if is_amplicon_read(read, track):
            call.n_contaminant += 1
            call.n_amplicon += 1
            continue
        ext = extract_barcode(read, track)
        if ext is not None:
            call.extractions.append(ext)
    return call


# --------------------------------------------------------------------------- #
# Bridge: called bc1 -> designed donor (magestic bc1->donor table)             #
# --------------------------------------------------------------------------- #
#
# Shengdi's R bridges with a pure exact ``left_join(by="bc1")`` (many-to-many,
# "take best hit"). Exact lookup is the workhorse; an optional HD1 neighbour probe
# and an in-set-HD1 collision flag are thin layers on top (a single sequencing
# error could have come from a different real barcode: random N20 => ~27 HD1
# collision pairs per 1e6 barcodes).

BRIDGE_EXACT = "EXACT"
BRIDGE_HD1 = "HD1"
BRIDGE_AMBIGUOUS = "AMBIGUOUS"
BRIDGE_NONE = "NONE"

_BASES = ("A", "C", "G", "T")


def hd1_neighbors(bc: str):
    """Yield the 3*len(bc) Hamming-distance-1 neighbours of ``bc``."""
    for i, ch in enumerate(bc):
        for b in _BASES:
            if b != ch:
                yield bc[:i] + b + bc[i + 1:]


@dataclass
class BridgeResult:
    """Outcome of bridging one called bc1 to the designed donor table."""

    barcode: str  # canonical query bc1
    tier: str  # EXACT / HD1 / AMBIGUOUS / NONE
    oligo: Optional[str]  # resolved best oligo (None when NONE / HD1-multi)
    support: int  # bc1->donor support of the resolved hit (0 if unknown)
    matched_barcode: Optional[str]  # table bc1 actually matched (== barcode for EXACT)
    n_oligos: int  # distinct oligos for the matched bc1 (source ambiguity)
    has_inset_hd1_neighbor: bool  # EXACT hit, but >=1 HD1 neighbour also in the table
    flags: Tuple[str, ...] = ()


class BridgeIndex:
    """bc1 -> designed-donor index over a magestic bc1->donor table.

    ``mapping`` is ``bc1 -> {oligo: summed_support}`` (a bc1 may appear against
    several oligos in the table; that is *source* ambiguity, kept and flagged).
    bc1 keys are stored in canonical orientation (== the caller's output); use
    ``table_revcomp`` in :meth:`from_tsv` if the table's bc1 column is reversed.
    """

    def __init__(self, mapping: Dict[str, Dict[str, int]]):
        self._map = mapping

    def __len__(self) -> int:
        return len(self._map)

    def __contains__(self, bc: str) -> bool:
        return bc in self._map

    @classmethod
    def from_records(cls, records) -> "BridgeIndex":
        """Build from an iterable of ``(bc1, oligo, support)`` (support may be 0)."""
        mapping: Dict[str, Dict[str, int]] = {}
        for bc1, oligo, support in records:
            bc1 = bc1.strip().upper()
            oligo = (oligo or "").strip()
            if not bc1 or not oligo:
                continue
            d = mapping.setdefault(bc1, {})
            d[oligo] = d.get(oligo, 0) + int(support or 0)
        return cls(mapping)

    @classmethod
    def from_tsv(
        cls,
        path: str,
        *,
        bc1_col: str = "bc1",
        oligo_col: str = "oligo_name",
        count_col: Optional[str] = None,
        min_count: int = 0,
        table_revcomp: bool = False,
    ) -> "BridgeIndex":
        """Load a bc1->donor table (header'd TSV). Aggregates support per (bc1, oligo).

        NB this loads the whole table into memory; for the multi-GB magestic tables
        run it on a compute node (the NNS bc1 reference is ~2 M rows / 3 GB). Uses an
        index-based ``csv.reader`` (not ``DictReader``) for speed/memory at that scale.
        """
        import csv

        mapping: Dict[str, Dict[str, int]] = {}
        with open(path, newline="") as fh:
            reader = csv.reader(fh, delimiter="\t")
            header = next(reader)
            idx = {name: i for i, name in enumerate(header)}
            bi, oi = idx[bc1_col], idx[oligo_col]
            ci = idx[count_col] if count_col else None
            for row in reader:
                if len(row) <= max(bi, oi):
                    continue
                bc = row[bi].strip().upper()
                oligo = row[oi].strip()
                if not bc or not oligo:
                    continue
                cnt = 0
                if ci is not None and ci < len(row):
                    cval = row[ci].strip()
                    if cval.isdigit():
                        cnt = int(cval)
                if cnt < min_count:
                    continue
                if table_revcomp:
                    bc = revcomp(bc)
                d = mapping.setdefault(bc, {})
                d[oligo] = d.get(oligo, 0) + cnt
        return cls(mapping)

    @staticmethod
    def _resolve(oligos: Dict[str, int]) -> Tuple[str, int]:
        """Best oligo by support; deterministic tie-break on the oligo name."""
        oligo, support = sorted(oligos.items(), key=lambda kv: (-kv[1], kv[0]))[0]
        return oligo, support

    def _has_inset_hd1(self, bc: str) -> bool:
        return any(nb in self._map for nb in hd1_neighbors(bc))

    def bridge(self, barcode: str, *, allow_hd1: bool = True) -> BridgeResult:
        bc = barcode.upper()
        if bc in self._map:
            oligos = self._map[bc]
            oligo, support = self._resolve(oligos)
            n = len(oligos)
            inset = self._has_inset_hd1(bc)
            flags = ("INSET_HD1_NEIGHBOR",) if inset else ()
            tier = BRIDGE_AMBIGUOUS if n > 1 else BRIDGE_EXACT
            return BridgeResult(bc, tier, oligo, support, bc, n, inset, flags)
        if allow_hd1:
            matches = {nb for nb in hd1_neighbors(bc) if nb in self._map}
            if len(matches) == 1:
                mb = next(iter(matches))
                oligos = self._map[mb]
                oligo, support = self._resolve(oligos)
                n = len(oligos)
                tier = BRIDGE_AMBIGUOUS if n > 1 else BRIDGE_HD1
                return BridgeResult(bc, tier, oligo, support, mb, n, False, ("HD1",))
            if len(matches) > 1:
                return BridgeResult(bc, BRIDGE_AMBIGUOUS, None, 0, None, len(matches), False, ("HD1_MULTI",))
        return BridgeResult(bc, BRIDGE_NONE, None, 0, None, 0, False, ())


# --------------------------------------------------------------------------- #
# Per-sample call: pick the top barcode, score confidence, bridge, emit        #
# --------------------------------------------------------------------------- #
#
# Confidence is feature-based (NOT a raw read count). On the 576-well validation
# the two-orientation agreement (both flank anchors extract the same 20mer) was a
# 100%-precision signal, and >=2 deduped molecules was 91-97%; a single-read
# single-anchor call was only 58%. So:
#   caller HIGH  <=> both_anchors OR support >= 2
# and the combined confidence folds in the bridge tier.


@dataclass
class TopBarcode:
    """The winning barcode of one track, with its confidence features."""

    barcode: str  # canonical
    support: int  # deduped molecules
    n_distinct: int  # distinct candidate barcodes seen on the track
    both_anchors: bool  # left- AND right-anchored reads agree (two-orientation)
    min_flank_mm: int  # cleanest anchored-flank mismatch count


def top_barcode(call: TrackCall) -> Optional[TopBarcode]:
    """Deterministic winner of a track: max deduped support, then fewest flank
    mismatches, then lexicographic barcode."""
    support = call.barcode_support()
    if not support:
        return None
    best_mm: Dict[str, int] = {}
    anchors: Dict[str, set] = {}
    for e in call.extractions:
        best_mm[e.barcode] = min(best_mm.get(e.barcode, 999), e.flank_mismatches)
        anchors.setdefault(e.barcode, set()).add(e.anchor)
    bc, sup = sorted(support.items(), key=lambda kv: (-kv[1], best_mm[kv[0]], kv[0]))[0]
    return TopBarcode(
        barcode=bc, support=sup, n_distinct=len(support),
        both_anchors={"left", "right"} <= anchors[bc], min_flank_mm=best_mm[bc],
    )


def caller_confidence(top: TopBarcode) -> str:
    """'high' if two-orientation agreement or >=2 deduped molecules, else 'low'."""
    return "high" if (top.both_anchors or top.support >= 2) else "low"


def _combine_confidence(caller: str, bridge_tier: str) -> str:
    """Fold caller confidence and bridge tier into one transparent label."""
    if bridge_tier in (BRIDGE_NONE, BRIDGE_AMBIGUOUS):
        return "low"
    if caller == "high" and bridge_tier == BRIDGE_EXACT:
        return "high"
    return "medium"


@dataclass
class Bc1Call:
    """The final per-sample bc1 call (caller features + donor bridge + confidence)."""

    sample: str
    track: str  # "genome" / "plasmid" / "none"
    barcode: str  # canonical bc1 ("" if no call)
    support: int
    n_distinct: int
    both_anchors: bool
    min_flank_mm: int
    n_reads_seen: int
    n_contaminant: int
    n_amplicon: int
    caller_confidence: str  # "high" / "low" / "none"
    bridge_tier: str
    oligo: Optional[str]
    bridge_support: int
    n_oligos: int
    bridge_flags: Tuple[str, ...]
    confidence: str  # combined high / medium / low / none

    @staticmethod
    def empty(sample: str) -> "Bc1Call":
        return Bc1Call(sample, "none", "", 0, 0, False, -1, 0, 0, 0, "none",
                       BRIDGE_NONE, None, 0, 0, (), "none")


def call_sample(
    bam: "pysam.AlignmentFile",
    sample: str,
    genome_track: TrackConfig,
    plasmid_track: Optional[TrackConfig] = None,
    bridge: Optional[BridgeIndex] = None,
) -> Bc1Call:
    """Call the bc1 for one sample BAM. Genome track is primary; the plasmid track
    is used only as a fallback **and only for a high-confidence** plasmid call
    (the plasmid track has no amplicon filter, so a low-confidence plasmid call is
    usually contamination — validated: naive genome-else-plasmid hurts precision)."""
    gcall = extract_track(bam, genome_track)
    top = top_barcode(gcall)
    track, used = "genome", gcall
    if top is None and plasmid_track is not None:
        pcall = extract_track(bam, plasmid_track)
        ptop = top_barcode(pcall)
        if ptop is not None and caller_confidence(ptop) == "high":
            top, track, used = ptop, "plasmid", pcall
    if top is None:
        return Bc1Call.empty(sample)

    conf = caller_confidence(top)
    br = bridge.bridge(top.barcode) if bridge is not None else BridgeResult(
        top.barcode, BRIDGE_NONE, None, 0, None, 0, False, ()
    )
    return Bc1Call(
        sample=sample, track=track, barcode=top.barcode, support=top.support,
        n_distinct=top.n_distinct, both_anchors=top.both_anchors, min_flank_mm=top.min_flank_mm,
        n_reads_seen=used.n_reads_seen, n_contaminant=used.n_contaminant, n_amplicon=used.n_amplicon,
        caller_confidence=conf, bridge_tier=br.tier, oligo=br.oligo, bridge_support=br.support,
        n_oligos=br.n_oligos, bridge_flags=br.flags, confidence=_combine_confidence(conf, br.tier),
    )


# --------------------------------------------------------------------------- #
# Config -> tracks, batch runner, CLI                                          #
# --------------------------------------------------------------------------- #


def load_locus_config(path: str) -> Dict[str, str]:
    """Parse a minimal ``key: 'value'`` config (the magestic ``input.yaml`` subset:
    barcode_locus / barcode_left_sequence / barcode_right_sequence and the plasmid_*
    equivalents). No PyYAML dependency."""
    cfg: Dict[str, str] = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.lstrip().startswith("#") or ":" not in line:
                continue
            key, val = line.split(":", 1)
            cfg[key.strip()] = val.strip().strip("'\"")
    return cfg


def _parse_locus(s: str) -> Tuple[str, int, int]:
    contig, rng = s.split(":")
    a, b = rng.split("-")
    return contig, int(a), int(b)


def _parse_window(s: Optional[str]) -> Optional[Tuple[int, int]]:
    if not s:
        return None
    a, b = s.split("-")
    return int(a), int(b)


def tracks_from_config(
    cfg: Dict[str, str],
    *,
    amplicon_left_window: Optional[Tuple[int, int]] = None,
    amplicon_right_window: Optional[Tuple[int, int]] = None,
) -> Tuple[TrackConfig, TrackConfig]:
    """Build the (genome, plasmid) tracks from a locus config. Genome bc1 is read
    out as revcomp(canonical) so genome ``revcomp=True``; the plasmid cassette is
    inverted so plasmid-forward is already canonical (``revcomp=False``)."""
    gchr, gs, ge = _parse_locus(cfg["barcode_locus"])
    pchr, ps, pe = _parse_locus(cfg["plasmid_bc_locus"])
    genome = TrackConfig(
        name="genome", contig=gchr, bc_start=gs, bc_end=ge,
        left_junction=cfg["barcode_left_sequence"], right_junction=cfg["barcode_right_sequence"],
        revcomp=True, amplicon_left_window=amplicon_left_window, amplicon_right_window=amplicon_right_window,
    )
    plasmid = TrackConfig(
        name="plasmid", contig=pchr, bc_start=ps, bc_end=pe,
        left_junction=cfg["plasmid_bc_left_sequence"], right_junction=cfg["plasmid_bc_right_sequence"],
        revcomp=False,
    )
    return genome, plasmid


EMIT_COLUMNS = [
    "sample", "track", "barcode", "support", "n_distinct", "both_anchors", "min_flank_mm",
    "n_reads_seen", "n_contaminant", "n_amplicon", "caller_confidence",
    "bridge_tier", "oligo", "bridge_support", "n_oligos", "bridge_flags", "confidence",
]


def call_to_row(c: Bc1Call) -> List[str]:
    return [
        c.sample, c.track, c.barcode, str(c.support), str(c.n_distinct), str(int(c.both_anchors)),
        str(c.min_flank_mm), str(c.n_reads_seen), str(c.n_contaminant), str(c.n_amplicon),
        c.caller_confidence, c.bridge_tier, c.oligo or "", str(c.bridge_support), str(c.n_oligos),
        ";".join(c.bridge_flags), c.confidence,
    ]


def _sample_from_bam(path: str, suffix: str) -> str:
    base = os.path.basename(path)
    return re.sub(re.escape(suffix) + r"$", "", base)


def run(
    bam_dir: str,
    config: str,
    fasta: str,
    out: str,
    *,
    bridge_table: Optional[str] = None,
    bridge_bc1_col: str = "bc1",
    bridge_oligo_col: str = "oligo_name",
    bridge_count_col: Optional[str] = None,
    amplicon_left_window: Optional[str] = None,
    amplicon_right_window: Optional[str] = None,
    use_plasmid: bool = False,
    bam_suffix: str = ".sort.dedup.recal2.bam",
) -> int:
    cfg = load_locus_config(config)
    genome, plasmid = tracks_from_config(
        cfg,
        amplicon_left_window=_parse_window(amplicon_left_window),
        amplicon_right_window=_parse_window(amplicon_right_window),
    )
    with pysam.FastaFile(fasta) as fa:
        genome.validate_against_reference(fa)
        plasmid.validate_against_reference(fa)
    bridge = None
    if bridge_table:
        sys.stderr.write(f"[bridge] loading {bridge_table} ...\n")
        bridge = BridgeIndex.from_tsv(
            bridge_table, bc1_col=bridge_bc1_col, oligo_col=bridge_oligo_col, count_col=bridge_count_col,
        )
        sys.stderr.write(f"[bridge] {len(bridge)} distinct bc1\n")
    bams = sorted(glob.glob(os.path.join(bam_dir, "*" + bam_suffix)))
    sys.stderr.write(f"[run] {len(bams)} bams\n")
    with open(out, "w") as o:
        o.write("\t".join(EMIT_COLUMNS) + "\n")
        for i, bampath in enumerate(bams, 1):
            sample = _sample_from_bam(bampath, bam_suffix)
            with pysam.AlignmentFile(bampath, "rb") as af:
                call = call_sample(af, sample, genome, plasmid if use_plasmid else None, bridge)
            o.write("\t".join(call_to_row(call)) + "\n")
            if i % 100 == 0:
                sys.stderr.write(f"  ...{i}/{len(bams)}\n")
    sys.stderr.write(f"[done] {out}\n")
    return 0


def _build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Call bc1 lineage barcodes from WGS BAMs and bridge to designed donors "
        "(TRACE WGS Stage 2-ii).",
    )
    p.add_argument("--bam-dir", required=True, help="directory of per-sample BAMs")
    p.add_argument("--config", required=True, help="locus config (input.yaml subset)")
    p.add_argument("--fasta", required=True, help="reference fasta with both barcode contigs (indexed)")
    p.add_argument("--out", required=True, help="output TSV of per-sample calls")
    p.add_argument("--bridge-table", help="bc1->donor TSV (e.g. the magestic bc1 reference table)")
    p.add_argument("--bridge-bc1-col", default="bc1")
    p.add_argument("--bridge-oligo-col", default="oligo_name")
    p.add_argument("--bridge-count-col", default=None)
    p.add_argument("--amplicon-left-window", default=None, help="e.g. 10224-10235 (enables amplicon filter)")
    p.add_argument("--amplicon-right-window", default=None, help="e.g. 10318-10326")
    p.add_argument("--use-plasmid", action="store_true", help="allow high-confidence plasmid-track fallback")
    p.add_argument("--bam-suffix", default=".sort.dedup.recal2.bam")
    return p


def main(argv: Optional[List[str]] = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    return run(
        args.bam_dir, args.config, args.fasta, args.out,
        bridge_table=args.bridge_table, bridge_bc1_col=args.bridge_bc1_col,
        bridge_oligo_col=args.bridge_oligo_col, bridge_count_col=args.bridge_count_col,
        amplicon_left_window=args.amplicon_left_window, amplicon_right_window=args.amplicon_right_window,
        use_plasmid=args.use_plasmid, bam_suffix=args.bam_suffix,
    )


if __name__ == "__main__":
    raise SystemExit(main())
