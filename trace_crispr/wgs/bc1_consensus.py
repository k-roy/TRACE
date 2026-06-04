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

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import pysam

from ..core.cigar import parse_cigar_to_operations

# Proper-pair primary flags Shengdi kept (read1/read2 x forward/reverse). Used as
# a feature, not a hard gate (we keep informative reads at n=1).
PROPER_PAIR_FLAGS = frozenset({83, 99, 147, 163})

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
        """Map a forward-on-contig barcode to the canonical bc1 orientation
        (genome-forward == the orientation of the bc1->donor reference table)."""
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
        if q is None:
            return None
        barcode = read.query_sequence[q : q + n]
    else:
        anchor = "right"
        q_end = ref_to_query_offset(read, bc_end0)  # first base past barcode
        if q_end is None:
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
    """Collapse PCR duplicates. Uses the fragment signature when the mate is
    mapped, else a single-read fallback (ref span + cigar) so true singletons
    are not inflated by duplicates."""
    if read.is_paired and not read.mate_is_unmapped and read.next_reference_id == read.reference_id:
        return ("pair", min(read.reference_start, read.next_reference_start), read.template_length, read.is_read1)
    return ("single", read.reference_start, read.reference_end, read.cigarstring, read.is_read1)


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
