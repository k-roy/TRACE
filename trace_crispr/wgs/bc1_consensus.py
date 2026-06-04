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


@dataclass
class TrackConfig:
    """One barcode track (genomic-integrated or episomal-plasmid)."""

    name: str  # e.g. "genome" / "plasmid"
    contig: str
    bc_start: int  # 1-based inclusive (first barcode base)
    bc_end: int  # 1-based inclusive (last barcode base)
    left_junction: str  # constant left flank (genome orientation)
    right_junction: str  # constant right flank

    @property
    def bc_len(self) -> int:
        return self.bc_end - self.bc_start + 1

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
    """A barcode extracted from one read, with raw per-read features."""

    barcode: str
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
    barcode = barcode.upper()

    flank_mm = _flank_mismatches(read, track, anchor)
    return Extraction(
        barcode=barcode,
        anchor=anchor,
        read_name=read.query_name,
        is_reverse=bool(read.is_reverse),
        mapq=read.mapping_quality,
        flank_mismatches=flank_mm,
        barcode_has_n=("N" in barcode),
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
    """Read-2 amplicon-hop contaminant test (Kevin's experimental signature):
    a barcode-bearing read whose mate maps **off** the barcode contig. For a
    genuine integration read both ends sit near the locus (mate on-contig); a
    hopped bc1-amplicon read-2 has its genuine WGS read-1 mate on a real
    chromosome. Drop on mate-off-contig alone (the agreed high-recall rule)."""
    if read.mate_is_unmapped:
        return True
    return read.next_reference_name != track.contig


@dataclass
class TrackCall:
    """All extractions for one (sample, track), pre-bridge."""

    track: str
    extractions: List[Extraction] = field(default_factory=list)
    n_reads_seen: int = 0
    n_contaminant: int = 0

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
    start = max(0, track.bc_start - 1 - window_pad)
    end = track.bc_end + window_pad
    for read in bam.fetch(track.contig, start, end):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        # only reads that actually overlap the barcode interval can carry it
        if read.reference_end <= track.bc_start - 1 or read.reference_start >= track.bc_end:
            continue
        call.n_reads_seen += 1
        # Primary contaminant defense: a genuine barcode read is a PROPER PAIR
        # (both ends near the locus). A read-2 amplicon hop is not — its genuine
        # WGS read-1 mate sits on a real chromosome — so requiring is_proper_pair
        # drops the hops (subsumes the mate-off-contig rule).
        if not read.is_proper_pair or is_read2_hop(read, track):
            call.n_contaminant += 1
            continue
        ext = extract_barcode(read, track)
        if ext is not None:
            call.extractions.append(ext)
    return call
