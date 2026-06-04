"""
trace_crispr.wgs.sv_plasmid
===========================

Structural-variant (SV) and plasmid-integration detection on a **combined-
reference** coordinate-sorted BAM of pooled MAGESTIC WGS data.

The combined reference is the genome **plus** per-oligo ``plasmid_{oligo}`` and
``integrated_plasmid_{oligo}`` contigs (see :func:`build_integrated_plasmid_contig`).
On that reference this module does three complementary things at each target
(donor) window:

* **Half B — donor-read filter + chimeric/translocation SV** (the heart):
  :func:`clean_donor_reads` reproduces Shengdi Li's ``clean_donor_reads.pl``
  bit-for-bit. It partitions the read pairs overlapping the donor window into
  *genome* (real target reads), *plasmid* (donor-template / intact-plasmid reads
  that must be discarded before variant calling) and *chimeric* (potential
  REDI-barcode-locus → target translocations), writes the per-sample
  ``*.target_filtering.anno`` fragment-count table, and emits a cleaned BAM with
  the plasmid/unresolved reads removed from the target site.

* **Half A — plasmid-integration detection** (Kevin R. Roy): :func:`detect_plasmid_integration`
  counts primary reads that span the genome↔integration junctions of each
  ``integrated_plasmid_{oligo}`` contig — direct evidence that the whole plasmid
  integrated at the cut site rather than being resolved as a clean edit.

* **Half C — mate-geometry large deletion** (Kevin R. Roy):
  :func:`one_read_inside_large_deletion` flags read pairs whose fragment span
  exceeds the nested mate by >1 kb (a large genomic deletion spanned by the pair);
  :func:`within_read_large_deletions` reuses :mod:`trace_crispr.core.cigar` for
  within-read large deletions.

Provenance
----------
Half B is a faithful Python port of **Shengdi Li**'s ``clean_donor_reads.pl``
(2022-10-30), in ``WGS_MAGESTIC_QTL_outcome_mapping`` — part of the MAGESTIC-SCORE
workflow (https://github.com/shli-embl/MAGESTIC-SCORE). Halves A/C and the
reference builder are ports of **Kevin R. Roy**'s
``screen_bam_for_plasmid_integration.py`` /
``construct_circular_and_integrated_plasmid_fasta.py`` /
``pysam_helper_functions.combine_aligned_segments_and_metrics``.

Parity fidelity (Half B)
------------------------
Validated bit-for-bit against Shengdi's reference ``*.target_filtering.anno``
outputs. The port reproduces the Perl semantics **exactly**, including several
deliberate quirks that matter for parity:

* read reference span is **M+D only** (``length_CIGAR``); ``N``/``I``/``S``/``H``
  do not extend it (see :func:`ref_span_md`), unlike pysam ``reference_end``;
* every alignment record is processed (the Perl only skips ``@`` header lines) —
  secondary/supplementary reads are **not** filtered out in Half B;
* the duplicate gate is Shengdi's buggy binary-``substr`` hack, not ``flag & 0x400``
  (see :func:`decode_duplicate`; ``mode="parity"`` reproduces it);
* the near-mate path tests *fragment* containment with **no** ``-1`` while the
  distal path tests *read* containment **with** ``-1`` — an internal coordinate
  inconsistency that is reproduced verbatim;
* the variant-midpoint coverage ``v_coord = (donor_start + donor_end) / 2`` is a
  Python ``float`` (half-integer when the sum is odd) compared with strict ``<``/``>``;
* a read pair may land in more than one fragment category at once (counted in
  every column it reaches) and *discard beats include* for cleaned-BAM membership.
"""

from __future__ import annotations

import argparse
import os
import re
import subprocess
from dataclasses import dataclass, field
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Set, Tuple, TypedDict, Union

from ..core.cigar import get_deletions_from_cigar

# Default flexible boundary (Shengdi's ``$flexible_bases``): the donor interval is
# widened by this many bp before containment tests, because a construct has a
# chance to match genome sequence; 6 bp is the maximum expected match by chance.
FLEXIBLE_BASES = 6

_CIGAR_SPAN_RE = re.compile(r"(\d+)([MIDSH])")
_LEADING_SOFTCLIP_RE = re.compile(r"^\d+S")
_TRAILING_SOFTCLIP_RE = re.compile(r"\d+S$")


# --------------------------------------------------------------------------- #
# CIGAR span + duplicate decoding (the two fidelity-critical primitives)       #
# --------------------------------------------------------------------------- #


def ref_span_md(cigar: str) -> int:
    """Reference span consumed by a read = sum of ``M`` and ``D`` lengths.

    This is a faithful port of Shengdi's ``length_CIGAR``: only ``M`` and ``D``
    extend the span; ``I``/``S``/``H`` do not, and the Perl regex does not even
    recognise ``N`` (skipped region) — it would die. WGS DNA alignments contain
    no ``N``, so we raise :class:`ValueError` to surface any unexpected ``N`` the
    same way the original would have failed, rather than silently miscounting
    (which is what pysam ``reference_end`` would do, since it counts ``N``).

    A bare ``*`` (unavailable CIGAR) is treated as span 0, as in the Perl.
    """
    if cigar == "*":
        return 0
    span = 0
    consumed = 0
    for length, op in _CIGAR_SPAN_RE.findall(cigar):
        consumed += len(length) + 1
        if op in ("M", "D"):
            span += int(length)
    if consumed != len(cigar):
        raise ValueError(f"CIGAR string unrecognized (only IDMSH supported): {cigar!r}")
    return span


def decode_duplicate(flag: int, mode: str = "standard") -> bool:
    """Decide whether a SAM ``FLAG`` marks a PCR/optical duplicate.

    ``mode="standard"`` (production default): the correct ``flag & 0x400`` test.

    ``mode="parity"``: reproduces Shengdi's ``if_duplicate`` **bit-for-bit**, which
    converts the flag to a big-endian binary string (leading zeros stripped) and
    reads ``substr(str, 10, 1)`` — the 11th character *from the left*. Because the
    string is right-aligned, that index does **not** address the 0x400 duplicate
    bit for most flags (e.g. ``1024`` → ``"10000000000"`` → index 10 is the
    least-significant ``paired`` bit), so the historical pipeline mislabels
    duplicates. The cleaned-BAM membership is unaffected (it never consults the
    duplicate flag), but the ``*.target_filtering.anno`` fragment **counts** are,
    so reproducing this exactly is required to match Shengdi's reference outputs.
    """
    if mode == "standard":
        return bool(flag & 0x400)
    if mode != "parity":
        raise ValueError(f"unknown duplicate mode {mode!r} (expected 'standard' or 'parity')")
    binstr = format(flag, "b")
    if len(binstr) < 11:
        return False
    return binstr[10] == "1"


# --------------------------------------------------------------------------- #
# Interval predicates (faithful to the Perl helpers)                           #
# --------------------------------------------------------------------------- #


def _is_overlapped(l1: int, r1: int, l2: int, r2: int) -> bool:
    """``True`` unless the closed intervals ``[l1,r1]`` and ``[l2,r2]`` are disjoint."""
    return not (r1 < l2 or l1 > r2)


def _is_contained(l1: int, r1: int, l2: int, r2: int) -> bool:
    """``True`` if interval ``[l2,r2]`` fully contains ``[l1,r1]`` (Perl ``is_contained``)."""
    return l2 <= l1 and r2 >= r1


def _is_softclipped_near_junction(
    interval_left: int,
    interval_right: int,
    read_map_start: int,
    cigar: str,
    flexible_boundary: int,
) -> bool:
    """Port of ``is_softclipped_near_junction``.

    A read soft-clipped within ``flexible_boundary`` bp of a donor boundary is a
    candidate REDI↔target translocation. A leading soft-clip is tested against the
    donor left boundary; a trailing soft-clip against the right boundary. Leading
    takes precedence (Perl ``if`` / ``elsif``).
    """
    if _LEADING_SOFTCLIP_RE.match(cigar):
        return abs(read_map_start - interval_left) <= flexible_boundary
    if _TRAILING_SOFTCLIP_RE.search(cigar):
        return abs(read_map_start + ref_span_md(cigar) - interval_right) <= flexible_boundary
    return False


# --------------------------------------------------------------------------- #
# REDI / cassette black-list                                                    #
# --------------------------------------------------------------------------- #

# (chrom, start, end) genomic+contig regions where REDI-cassette reads mis-map.
# The blacklist is a REQUIRED co-equal input: it is a *mix* of genomic
# mis-mapping loci (marker/promoter regions on real chromosomes) and plasmid /
# cassette contigs, so the combined reference does not subsume it.
Blacklist = Dict[str, List[Tuple[int, int]]]


def _coerce_int(token: str) -> int:
    """Numeric coercion matching Perl's behaviour (``undef``/non-numeric → 0).

    The original parses the blacklist with ``split("\\t")`` and indexes columns
    0/1/2; a malformed row (e.g. the file's space-delimited
    ``yT177_after_bc1_integration`` line, which has no tabs) leaves start/end as
    ``undef``, which Perl coerces to ``0`` in the numeric ``start < coord < end``
    test — making that entry inert (it can never match). We reproduce that.
    """
    try:
        return int(token)
    except ValueError:
        return 0


def load_blacklist(path: str) -> Blacklist:
    """Load ``REDI_reads_mapped_regions.txt`` into ``{chrom: [(start, end), ...]}``.

    Rows are ``chrom<TAB>start<TAB>end`` with an optional trailing label column
    (ignored, as in the Perl). Parsing is intentionally permissive to match
    Shengdi's ``split("\\t")`` exactly, including the file's one inert
    space-delimited row (see :func:`_coerce_int`).
    """
    bl: Blacklist = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            cols = line.split("\t")
            chrom = cols[0]
            start = _coerce_int(cols[1]) if len(cols) > 1 else 0
            end = _coerce_int(cols[2]) if len(cols) > 2 else 0
            bl.setdefault(chrom, []).append((start, end))
    return bl


def _is_mapped_to_cassette(chrom: str, coord: int, blacklist: Blacklist) -> bool:
    """Port of ``is_mapped_to_cassete``: ``start < coord < end`` (strict both sides)."""
    for start, end in blacklist.get(chrom, ()):
        if start < coord < end:
            return True
    return False


# --------------------------------------------------------------------------- #
# Half B: donor-read filter + SV classification                                #
# --------------------------------------------------------------------------- #


@dataclass
class SamRecord:
    """The SAM columns Half B needs, parsed exactly as the Perl ``split("\\t")``.

    ``pos``/``pnext`` are 1-based (SAM ``POS``/``PNEXT`` verbatim); ``rnext`` keeps
    the literal ``=`` used for a same-reference mate.
    """

    qname: str
    flag: int
    rname: str
    pos: int
    cigar: str
    rnext: str
    pnext: int
    tlen: int


def parse_sam_line(line: str) -> SamRecord:
    """Parse one alignment line of ``samtools view`` output into a :class:`SamRecord`."""
    c = line.rstrip("\n").split("\t")
    return SamRecord(
        qname=c[0],
        flag=int(c[1]),
        rname=c[2],
        pos=int(c[3]),
        cigar=c[5],
        rnext=c[6],
        pnext=int(c[7]),
        tlen=int(c[8]),
    )


@dataclass
class CleanResult:
    """Outcome of classifying a sample's donor-overlapping reads (Half B pass 1)."""

    overlapped: Set[str] = field(default_factory=set)
    discarded: Set[str] = field(default_factory=set)
    included: Set[str] = field(default_factory=set)
    redi_translocation: Set[str] = field(default_factory=set)
    flag: Dict[str, int] = field(default_factory=dict)
    variant_cov: int = 0
    unique_genome: Set[str] = field(default_factory=set)
    unique_plasmid: Set[str] = field(default_factory=set)
    unique_chimeric: Set[str] = field(default_factory=set)

    @property
    def frags_genome(self) -> int:
        return len(self.unique_genome)

    @property
    def frags_plasmid(self) -> int:
        return len(self.unique_plasmid)

    @property
    def frags_chimeric(self) -> int:
        return len(self.unique_chimeric)

    def anno_row(self) -> str:
        """The single data row of ``*.target_filtering.anno`` (col1 is variant_cov)."""
        return f"{self.variant_cov}\t{self.frags_genome}\t{self.frags_plasmid}\t{self.frags_chimeric}"


# The ``*.target_filtering.anno`` header is mislabeled in the original (column 1 is
# actually the variant-midpoint read count, not "read_count_mid_point" as a frag
# count); it is preserved verbatim because downstream code parses it positionally.
ANNO_HEADER = "read_count_mid_point\tfrags_genome\tfrags_plasmid\tfrags_chimeric"


def classify_donor_reads(
    records: Iterable[SamRecord],
    chrom: str,
    donor_start: int,
    donor_end: int,
    blacklist: Blacklist,
    flexible_bases: int = FLEXIBLE_BASES,
    dup_mode: str = "standard",
) -> CleanResult:
    """Classify every alignment record overlapping a donor window (Half B pass 1).

    Faithful port of the counting pass of ``clean_donor_reads.pl``. ``records``
    must yield the sample's alignment records (header lines already stripped).

    Args:
        records: iterable of :class:`SamRecord` (one per alignment line).
        chrom, donor_start, donor_end: target donor window (1-based inclusive).
        blacklist: REDI/cassette regions from :func:`load_blacklist`.
        flexible_bases: donor-boundary slack (Shengdi default 6).
        dup_mode: ``"parity"`` to reproduce Shengdi's fragment counts exactly,
            ``"standard"`` for the correct duplicate gate.

    Returns:
        A :class:`CleanResult` with the per-read-name fragment sets and counts.
    """
    res = CleanResult()
    v_coord = (donor_start + donor_end) / 2  # float midpoint, as in Perl

    for rec in records:
        span = ref_span_md(rec.cigar)
        read_end = rec.pos + span  # NB: no -1, matching the Perl overlap test

        if not (rec.rname == chrom and _is_overlapped(donor_start, donor_end, rec.pos, read_end)):
            continue

        res.overlapped.add(rec.qname)
        res.flag[rec.qname] = rec.flag
        is_dup = decode_duplicate(rec.flag, dup_mode)

        if rec.rnext == "=" and abs(rec.tlen) < 1000:
            # Mate maps near: judge by fragment span.
            frag_start = min(rec.pos, rec.pnext)
            frag_end = frag_start + abs(rec.tlen) - 1
            if _is_contained(frag_start, frag_end, donor_start - flexible_bases, donor_end + flexible_bases):
                # Fragment lies within the donor → plasmid/donor-template read; discard.
                res.discarded.add(rec.qname)
                if not is_dup:
                    res.unique_plasmid.add(rec.qname)
            else:
                res.included.add(rec.qname)
                if _is_softclipped_near_junction(donor_start, donor_end, rec.pos, rec.cigar, flexible_bases):
                    res.redi_translocation.add(rec.qname)
                    if not is_dup:
                        res.unique_chimeric.add(rec.qname)
                else:
                    if not is_dup:
                        res.unique_genome.add(rec.qname)
                if not is_dup and rec.pos < v_coord < read_end:
                    res.variant_cov += 1

        elif rec.rnext != "*":
            # Mate maps distal (or to another chromosome).
            mate_chr = rec.rname if rec.rnext == "=" else rec.rnext
            if _is_mapped_to_cassette(mate_chr, rec.pnext, blacklist):
                # NB: distal path tests *read* containment WITH -1 (Perl quirk).
                if not _is_contained(
                    rec.pos, rec.pos + span - 1, donor_start - flexible_bases, donor_end + flexible_bases
                ):
                    res.included.add(rec.qname)
                    res.redi_translocation.add(rec.qname)
                    if not is_dup:
                        res.unique_chimeric.add(rec.qname)
                else:
                    res.discarded.add(rec.qname)
                    if not is_dup:
                        res.unique_plasmid.add(rec.qname)
            else:
                res.included.add(rec.qname)
                if not is_dup:
                    res.unique_genome.add(rec.qname)
                    if rec.pos < v_coord < read_end:
                        res.variant_cov += 1
        else:
            # Mate not properly mapped → discard.
            res.discarded.add(rec.qname)

    return res


def _iter_sam(bam_path: str, samtools: str = "samtools") -> Iterator[str]:
    """Yield SAM lines (incl. ``@`` header) from ``samtools view -h <bam>``.

    Raises :class:`subprocess.CalledProcessError` if ``samtools`` exits non-zero
    (e.g. a truncated or unreadable BAM), so callers never silently process a
    partial SAM stream and emit wrong counts. The exit code is only checked on
    normal exhaustion of the generator; a consumer that stops early triggers
    ``GeneratorExit`` and skips the check (no spurious error on a broken pipe).
    """
    proc = subprocess.Popen(
        [samtools, "view", "-h", bam_path],
        stdout=subprocess.PIPE,
        text=True,
    )
    if proc.stdout is None:  # pragma: no cover - PIPE always provides a stdout
        raise RuntimeError("failed to capture samtools stdout")
    try:
        for line in proc.stdout:
            yield line
    finally:
        proc.stdout.close()
        rc = proc.wait()
    if rc:
        raise subprocess.CalledProcessError(rc, proc.args)


def write_anno(result: CleanResult, path: str) -> None:
    """Write the two-line ``*.target_filtering.anno`` table (header verbatim)."""
    with open(path, "w") as oh:
        oh.write(ANNO_HEADER + "\n")
        oh.write(result.anno_row() + "\n")


def clean_donor_reads(
    in_bam: str,
    out_bam: str,
    chrom: str,
    donor_start: int,
    donor_end: int,
    blacklist_path: str,
    flexible_bases: int = FLEXIBLE_BASES,
    dup_mode: str = "standard",
    samtools: str = "samtools",
) -> CleanResult:
    """Filter donor/plasmid/unresolved reads from a target site (Half B driver).

    Two passes over the input BAM (exactly like the Perl):

    1. classify every donor-overlapping read pair and write
       ``{out_bam}.target_filtering.anno``;
    2. write ``out_bam`` keeping all non-target reads plus the donor-overlapping
       reads that were *included* and not *discarded* (**discard beats include**).

    Cleaned-BAM membership is independent of ``dup_mode``; only the ``.anno``
    fragment counts depend on it (use ``dup_mode="parity"`` to reproduce Shengdi's
    reference counts).

    Returns the :class:`CleanResult` from pass 1.
    """
    blacklist = load_blacklist(blacklist_path)

    records = (parse_sam_line(ln) for ln in _iter_sam(in_bam, samtools) if not ln.startswith("@"))
    result = classify_donor_reads(
        records, chrom, donor_start, donor_end, blacklist, flexible_bases, dup_mode
    )
    write_anno(result, out_bam + ".target_filtering.anno")

    tmp_sam = out_bam + ".tmp.sam"
    with open(tmp_sam, "w") as out:
        for ln in _iter_sam(in_bam, samtools):
            if ln.startswith("@"):
                out.write(ln)
                continue
            name = ln.split("\t", 1)[0]
            if name in result.overlapped:
                if name in result.discarded:
                    continue
                if name in result.included:
                    out.write(ln)
                # else: unknown / unresolved → drop
            else:
                out.write(ln)
    with open(out_bam, "wb") as out_fh:
        subprocess.run([samtools, "view", "-bS", tmp_sam], stdout=out_fh, check=True)
    os.remove(tmp_sam)
    return result


# --------------------------------------------------------------------------- #
# Half A: plasmid-integration detection                                        #
# --------------------------------------------------------------------------- #


def _integrated_contigs(references: Sequence[str]) -> List[str]:
    return [r for r in references if r.startswith("integrated_plasmid_")]


def detect_plasmid_integration(
    bam_path: str,
    integrated_contigs: Optional[Sequence[str]] = None,
    junction_distance: Union[int, Dict[str, int]] = 1129,
    window: int = 200,
) -> Dict[str, int]:
    """Count primary reads spanning the genome↔integration junctions of a sample.

    Faithful port of Kevin R. Roy's ``screen_bam_for_plasmid_integration.py``. For
    each ``integrated_plasmid_{oligo}`` contig, primary reads (not unmapped,
    secondary or supplementary; **no** dedup, **no** ``NM`` filter) whose alignment
    spans the upstream junction (at ``junction_distance``) or the downstream
    junction (at ``contig_length - junction_distance``) within a ``±window`` fetch
    window are counted as plasmid-integration evidence.

    ``junction_distance`` is the parity literal (1129); for production it should be
    derived per contig as ``upstream_arm_length + donor_length`` (see
    :func:`integrated_plasmid_junction`). Pass a ``{contig: distance}`` mapping to
    set it per contig.

    Returns ``{"upstream_count", "downstream_count"}`` summed over all contigs for
    this (single-sample) BAM.
    """
    import pysam

    bam = pysam.AlignmentFile(bam_path, "rb")
    try:
        refs = bam.references
        contigs = list(integrated_contigs) if integrated_contigs is not None else _integrated_contigs(refs)
        up = 0
        down = 0
        for ref in contigs:
            if ref not in refs:
                continue
            if isinstance(junction_distance, dict):
                jd = junction_distance.get(ref)
                if jd is None:
                    continue
            else:
                jd = junction_distance
            up_junction = jd
            down_junction = bam.get_reference_length(ref) - jd
            for junction, side in ((up_junction, "up"), (down_junction, "down")):
                start = max(0, junction - window)
                end = junction + window
                for read in bam.fetch(ref, start, end):
                    if read.is_unmapped or read.is_secondary or read.is_supplementary:
                        continue
                    if read.reference_start <= junction < read.reference_end:
                        if side == "up":
                            up += 1
                        else:
                            down += 1
        return {"upstream_count": up, "downstream_count": down}
    finally:
        bam.close()


class PlasmidIntegrationRow(TypedDict):
    """One row of the plasmid-integration summary (the ``plasmid_integration_summary.tsv`` schema)."""

    sample_name: str
    upstream_count: int
    downstream_count: int
    total_count: int


def summarize_plasmid_integration(
    sample_bams: Dict[str, str],
    integrated_contigs: Optional[Sequence[str]] = None,
    junction_distance: Union[int, Dict[str, int]] = 1129,
    window: int = 200,
) -> List[PlasmidIntegrationRow]:
    """Run :func:`detect_plasmid_integration` over ``{sample_name: bam_path}``.

    Returns one :class:`PlasmidIntegrationRow` per sample, sorted by ``total_count``
    descending — the row schema of Kevin's ``plasmid_integration_summary.tsv``.
    """
    rows: List[PlasmidIntegrationRow] = []
    for sample, bam_path in sample_bams.items():
        counts = detect_plasmid_integration(bam_path, integrated_contigs, junction_distance, window)
        total = counts["upstream_count"] + counts["downstream_count"]
        rows.append(
            {
                "sample_name": sample,
                "upstream_count": counts["upstream_count"],
                "downstream_count": counts["downstream_count"],
                "total_count": total,
            }
        )
    rows.sort(key=lambda r: r["total_count"], reverse=True)
    return rows


# --------------------------------------------------------------------------- #
# Half C: mate-geometry + within-read large deletions                          #
# --------------------------------------------------------------------------- #


def one_read_inside_large_deletion(
    r1_start: int,
    r1_end: int,
    r2_start: int,
    r2_end: int,
    min_gap: int = 1000,
) -> bool:
    """Mate-geometry test for a large deletion spanned by a read pair.

    ``True`` when one mate is strictly nested inside the other **and** the pair's
    reference span exceeds the nested mate's span by more than ``min_gap`` bp — the
    geometry of a >1 kb genomic deletion bridged by the fragment. Coordinates are
    aligned reference positions (0-based, first/last aligned base of each read).

    Port of the ``one_read_inside_large_deletion`` branch of Kevin R. Roy's
    ``pysam_helper_functions.combine_aligned_segments_and_metrics``.
    """
    region_length = max(r1_end, r2_end) - min(r1_start, r2_start) + 1
    r1_length = r1_end - r1_start + 1
    r2_length = r2_end - r2_start + 1
    if region_length > r2_length + min_gap and (r2_start > r1_start and r2_end < r1_end):
        return True
    if region_length > r1_length + min_gap and (r1_start > r2_start and r1_end < r2_end):
        return True
    return False


def mate_pair_large_deletion(read_1, read_2, min_gap: int = 1000) -> bool:
    """pysam wrapper for :func:`one_read_inside_large_deletion`.

    Uses each read's first/last aligned reference position
    (``get_reference_positions()[0]`` / ``[-1]``), matching the original.
    """
    p1 = read_1.get_reference_positions()
    p2 = read_2.get_reference_positions()
    if not p1 or not p2:
        return False
    return one_read_inside_large_deletion(p1[0], p1[-1], p2[0], p2[-1], min_gap)


def within_read_large_deletions(read, cutoff: int = 50):
    """Within-read deletions strictly larger than ``cutoff`` bp.

    Reuses :func:`trace_crispr.core.cigar.get_deletions_from_cigar`. The original
    classifies a deletion as "large" when ``len(deleted_bases) > cutoff`` (strict),
    i.e. size ``>= cutoff + 1``.
    """
    return get_deletions_from_cigar(read, min_size=cutoff + 1)


# --------------------------------------------------------------------------- #
# Combined-reference construction                                              #
# --------------------------------------------------------------------------- #

_REVCOMP_TABLE = str.maketrans(
    "ACGTacgtRYMKrymkVBHDvbhdNn ", "TGCAtgcaYRKMyrkmBVDHbvdhNn "
)


def rev_comp(seq: str) -> str:
    """Reverse complement, preserving IUPAC ambiguity codes (Kevin's ``rev_comp``)."""
    return seq.translate(_REVCOMP_TABLE)[::-1]


def build_plasmid_contig(donor: str, backbone: str, guide: str, insert: str) -> str:
    """Linearised episomal plasmid contig ``donor + backbone + guide + insert``."""
    return f"{donor}{backbone}{guide}{insert}"


def build_shifted_plasmid_contig(plasmid_sequence: str, backbone_length: int) -> str:
    """Rotate a circular plasmid by half its backbone so the backbone-internal
    junction is mappable (the ``plasmid_sequence_shifted`` of the original)."""
    half = backbone_length // 2
    return plasmid_sequence[half:] + plasmid_sequence[:half]


def build_integrated_plasmid_contig(
    extended_donor: str,
    donor: str,
    integrated_segment: str,
    donor_strand: str,
    bp_to_extend: int = 1000,
) -> str:
    """Build an ``integrated_plasmid_{oligo}`` reference contig.

    Layout (genome → integrated plasmid → genome)::

        upstream_flank | donor | integrated_segment | donor | downstream_flank

    where the flanks are the outer ``bp_to_extend`` of ``extended_donor`` and
    ``integrated_segment = backbone + guide + insert`` (reverse-complemented when
    the donor is on the minus strand). Port of
    ``construct_integrated_plasmid_sequence``.
    """
    upstream = extended_donor[:bp_to_extend]
    downstream = extended_donor[-bp_to_extend:]
    segment = rev_comp(integrated_segment) if donor_strand == "-" else integrated_segment
    return f"{upstream}{donor}{segment}{donor}{downstream}"


def integrated_plasmid_junction(donor: str, bp_to_extend: int = 1000) -> int:
    """0-based contig coordinate of the upstream genome↔integration junction.

    Equals ``bp_to_extend + len(donor)`` (= ``junction_distance``); the downstream
    junction is symmetric at ``contig_length - junction``. For the parity library
    (``bp_to_extend=1000``, 129 bp donor) this is 1129.
    """
    return bp_to_extend + len(donor)


# --------------------------------------------------------------------------- #
# CLI (mirrors clean_donor_reads.pl)                                           #
# --------------------------------------------------------------------------- #


def _build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Filter donor-template / plasmid / chimeric reads from a target "
        "site and write the *.target_filtering.anno table (Python port of Shengdi "
        "Li's clean_donor_reads.pl)."
    )
    p.add_argument("--in", dest="in_bam", required=True)
    p.add_argument("--out", dest="out_bam", required=True)
    p.add_argument("--chr", dest="chrom", required=True)
    p.add_argument("--d-start", dest="donor_start", type=int, required=True)
    p.add_argument("--d-end", dest="donor_end", type=int, required=True)
    p.add_argument("--blacklist", required=True, help="REDI_reads_mapped_regions.txt")
    p.add_argument("--flexible-bases", type=int, default=FLEXIBLE_BASES)
    p.add_argument(
        "--dup-mode",
        choices=("standard", "parity"),
        default="standard",
        help="'parity' reproduces Shengdi's fragment counts bit-for-bit.",
    )
    p.add_argument("--samtools", default="samtools")
    return p


def main(argv: Optional[List[str]] = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    clean_donor_reads(
        in_bam=args.in_bam,
        out_bam=args.out_bam,
        chrom=args.chrom,
        donor_start=args.donor_start,
        donor_end=args.donor_end,
        blacklist_path=args.blacklist,
        flexible_bases=args.flexible_bases,
        dup_mode=args.dup_mode,
        samtools=args.samtools,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
