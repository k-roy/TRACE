"""
trace_crispr.wgs.sv_plasmid
===========================

Structural-variant (SV) and plasmid-integration detection on a **combined-
reference** coordinate-sorted BAM of single-sample MAGESTIC WGS data.

The combined reference is the genome **plus** per-oligo ``plasmid_{oligo}`` and
``integrated_plasmid_{oligo}`` contigs (see :func:`build_integrated_plasmid_contig`).
On that reference this module does three complementary things at each target
(donor) window:

* **Half B — donor-read filter + chimeric/translocation evidence**:
  :func:`clean_donor_reads` derives from Shengdi Li's ``clean_donor_reads.pl``.
  Standard mode corrects coordinate, filtering and category ambiguities; explicit
  ``dup_mode="parity"`` retains historical counting semantics. It partitions
  read pairs overlapping the donor window into
  *genome* (real target reads), *plasmid* (donor-template / intact-plasmid reads
  that must be discarded before variant calling) and *chimeric* (potential
  REDI-barcode-locus → target translocations), writes the per-sample
  ``*.target_filtering.anno`` fragment-count table, and emits a cleaned BAM with
  the plasmid/unresolved reads removed from the target site.

* **Half A — plasmid-integration detection** (Kevin R. Roy): :func:`detect_plasmid_integration`
  counts primary reads that span the genome↔integration junctions of each
  ``integrated_plasmid_{oligo}`` contig. These are candidate integration-support
  counts; junction reads alone do not establish complete plasmid integration.

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

Historical parity (Half B, explicit ``dup_mode="parity"`` only)
------------------------
The parity path preserves the Perl semantics below for replication. Production
standard mode uses reference-validated blacklist coordinates, primary/QC/MAPQ
filters, aligned-base overlap and mutually exclusive fragment categories instead.

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
import hashlib
import os
import re
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field
from contextlib import contextmanager
from typing import Dict, Hashable, Iterable, Iterator, List, NamedTuple, Optional, Sequence, Set, Tuple, TypedDict, Union

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


def load_blacklist(path: str, mode: str = "standard") -> Blacklist:
    """Load ``REDI_reads_mapped_regions.txt`` into ``{chrom: [(start, end), ...]}``.

    Standard mode requires tab-separated, 1-based inclusive intervals. Historical
    permissive Perl parsing is available only with ``mode="parity"``.
    """
    if mode not in ("standard", "parity"):
        raise ValueError(f"unknown blacklist mode {mode!r}")
    bl: Blacklist = {}
    with open(path) as fh:
        for line_number, line in enumerate(fh, 1):
            line = line.rstrip("\n")
            if not line:
                continue
            cols = line.split("\t")
            if mode == "standard":
                if len(cols) < 3 or not cols[0].strip():
                    raise ValueError(f"{path}:{line_number}: expected chrom<TAB>start<TAB>end")
                try:
                    start, end = int(cols[1]), int(cols[2])
                except ValueError as exc:
                    raise ValueError(f"{path}:{line_number}: non-integer blacklist interval") from exc
                if start < 1 or end < start:
                    raise ValueError(f"{path}:{line_number}: invalid 1-based blacklist interval")
                bl.setdefault(cols[0], []).append((start, end))
                continue
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
    the literal ``=`` used for a same-reference mate. ``rg`` is the ``RG:Z`` read
    group, or ``""`` when the record declares none (Stage4's missing-RG namespace).
    """

    qname: str
    flag: int
    rname: str
    pos: int
    cigar: str
    rnext: str
    pnext: int
    tlen: int
    mapq: int = 60
    rg: str = ""
    mate_cigar: str = ""  # MC:Z, the mate's CIGAR when the aligner wrote one


def parse_sam_line(line: str) -> SamRecord:
    """Parse one alignment line of ``samtools view`` output into a :class:`SamRecord`.

    The read group is part of fragment identity, so a record with more than one
    ``RG`` tag, a non-``Z`` type or an empty ID is rejected rather than resolved
    by picking one. This applies in parity mode too: such a line is invalid SAM.
    """
    c = line.rstrip("\n").split("\t")
    if len(c) < 11:
        raise ValueError("SAM alignment must have at least 11 tab-separated columns")
    tags = [f for f in c[11:] if f.startswith("RG:")]
    rg = ""
    if len(tags) > 1:
        raise ValueError(f"alignment {c[0]!r} declares {len(tags)} RG read-group tags; exactly one or none is allowed")
    if tags:
        parts = tags[0].split(":", 2)
        if len(parts) != 3 or parts[1] != "Z" or not parts[2]:
            raise ValueError(f"alignment {c[0]!r} has a malformed RG read-group tag {tags[0]!r}; expected RG:Z:<id>")
        rg = parts[2]
    mc = [f[5:] for f in c[11:] if f.startswith("MC:Z:")]
    return SamRecord(
        qname=c[0],
        flag=int(c[1]),
        rname=c[2],
        pos=int(c[3]),
        cigar=c[5],
        rnext=c[6],
        pnext=int(c[7]),
        tlen=int(c[8]),
        mapq=int(c[4]),
        rg=rg,
        mate_cigar=mc[0] if len(mc) == 1 else "",
    )


class FragmentKey(NamedTuple):
    """Standard-mode Half B fragment identity, the same namespace Stage4 counts.

    A paired fragment is ``(read_group, qname)`` with ``record=None``, so both
    mates and their supplementary records share one disposition; mate number is
    deliberately not part of it. An unpaired alignment is its own observation, so
    ``record`` holds its ``(rname, pos, cigar, reverse)``: a rejected single read
    cannot remove a valid same-name pair, or an unrelated same-name record.
    """

    read_group: str
    qname: str
    record: Optional[Tuple[str, int, str, bool]] = None


def fragment_key(rec: SamRecord) -> FragmentKey:
    """The :class:`FragmentKey` of one record (standard mode, both BAM passes)."""
    if rec.flag & 0x1:
        return FragmentKey(rec.rg, rec.qname)
    return FragmentKey(rec.rg, rec.qname, (rec.rname, rec.pos, rec.cigar, bool(rec.flag & 0x10)))


@dataclass
class CleanResult:
    """Outcome of classifying a sample's donor-overlapping reads (Half B pass 1).

    Set members identify fragments: bare query names in explicit parity mode
    (the historical QNAME-only namespace) and :class:`FragmentKey` in standard mode.
    """

    overlapped: Set[Hashable] = field(default_factory=set)
    discarded: Set[Hashable] = field(default_factory=set)
    included: Set[Hashable] = field(default_factory=set)
    redi_translocation: Set[Hashable] = field(default_factory=set)
    flag: Dict[Hashable, int] = field(default_factory=dict)
    variant_cov: int = 0
    unique_genome: Set[Hashable] = field(default_factory=set)
    unique_plasmid: Set[Hashable] = field(default_factory=set)
    unique_chimeric: Set[Hashable] = field(default_factory=set)
    # Donor-contained fragments RETAINED because the caller declared the plasmids purged
    # (``keep_contained=True``). Kept apart from ``unique_genome`` so a genome count never
    # silently includes them, and absent from the historical ``.anno`` row.
    contained_kept: Set[Hashable] = field(default_factory=set)
    unique_contained_kept: Set[Hashable] = field(default_factory=set)

    @property
    def frags_genome(self) -> int:
        return len(self.unique_genome)

    @property
    def frags_contained_kept(self) -> int:
        return len(self.unique_contained_kept)

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


def _classify_donor_reads_parity(
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


def _cigar_blocks(pos: int, cigar: str) -> tuple[list[tuple[int, int]], int, bool, bool]:
    """Aligned 1-based inclusive blocks, inclusive end, terminal soft clips.

    D/N advance the reference but never contribute aligned-base coverage. =/X
    have the same coordinate semantics as M. H/P may surround terminal S.
    """
    if cigar == "*":
        return [], pos - 1, False, False
    parts = re.findall(r"(\d+)([MIDNSHP=X])", cigar)
    if not parts or "".join(n + op for n, op in parts) != cigar or any(int(n) <= 0 for n, _ in parts):
        raise ValueError(f"invalid CIGAR {cigar!r}")
    blocks: list[tuple[int, int]] = []
    cursor = pos
    for n, op in parts:
        size = int(n)
        if op in "M=X":
            blocks.append((cursor, cursor + size - 1))
        if op in "MDN=X":
            cursor += size
    terminal = [op for _, op in parts if op not in "HP"]
    return blocks, cursor - 1, terminal[0] == "S", terminal[-1] == "S"


def _primary_sam(rec: SamRecord, min_mapq: int) -> bool:
    return not (rec.flag & (0x4 | 0x100 | 0x200 | 0x800)) and min_mapq <= rec.mapq < 255


def classify_donor_reads(
    records: Iterable[SamRecord], chrom: str, donor_start: int, donor_end: int,
    blacklist: Blacklist, flexible_bases: int = FLEXIBLE_BASES,
    dup_mode: str = "standard", min_mapq: int = 20, *, keep_contained: bool = False,
    plasmid_contigs: Iterable[str] = (),
) -> CleanResult:
    """Classify donor-overlapping alignments, with explicit historical parity.

    Standard mode uses aligned-base overlap, primary/QC/MAPQ gates and valid
    mapped-mate coordinates. Fragment categories are exclusive: donor-template
    discard takes precedence over chimeric, then genome. Midpoint coverage counts
    qualifying retained genome alignment records, not D/N spans or molecules.
    MAPQ 20 is a conservative policy default, not a calibrated accuracy threshold.
    Standard fragments are :class:`FragmentKey` values, the (RG, QNAME) pairing
    Stage4 uses, so a disposition in one read group never reaches another.

    ``keep_contained=True`` is for a sample whose plasmids are known to be purged
    (see :func:`plasmid_load`). A same-chromosome fragment lying wholly inside the
    donor +/- ``flexible_bases`` can then only be genomic, so it is RETAINED (in the
    cleaned BAM, counted in ``contained_kept``) instead of discarded as a possible
    plasmid donor template. Two cases stay discarded even then: a contained fragment
    soft-clipped at a donor edge (the clip is a junction to something else) and a
    read whose mate is on a cassette contig. Contained fragments never enter
    ``unique_genome`` or ``variant_cov``. Standard mode only.

    ``plasmid_contigs`` (opt-in, default none) names plasmid contigs whose mates count like cassette mates: a
    donor read with its mate anywhere on one is plasmid template when contained and chimeric otherwise. Without
    it, a mate on a plasmid outside the REDI blacklist intervals falls through to genome (planning/328t found
    no such read in the NNS cohort, so the default is left unchanged).
    """
    if (any(type(value) is not int for value in (donor_start, donor_end, flexible_bases, min_mapq))
            or donor_start < 1 or donor_end < donor_start or flexible_bases < 0 or not 0 <= min_mapq <= 254):
        raise ValueError("invalid donor interval, flexible bases or MAPQ threshold")
    if type(keep_contained) is not bool:
        raise ValueError("keep_contained must be a bool")
    if keep_contained and dup_mode != "standard":
        raise ValueError("keep_contained is a standard-mode policy; parity mode reproduces the historical discard")
    plasmid_contigs = frozenset(plasmid_contigs)
    if plasmid_contigs and dup_mode != "standard":
        raise ValueError("plasmid_contigs is a standard-mode policy")
    if dup_mode == "parity":
        return _classify_donor_reads_parity(records, chrom, donor_start, donor_end, blacklist, flexible_bases, dup_mode)
    if dup_mode != "standard":
        raise ValueError(f"unknown duplicate mode {dup_mode!r}")
    res = CleanResult()
    midpoint = (donor_start + donor_end) / 2
    coverage: Dict[FragmentKey, int] = {}
    for rec in records:
        if rec.rname != chrom or rec.pos < 1:
            continue
        blocks, read_end, leading_soft, trailing_soft = _cigar_blocks(rec.pos, rec.cigar)
        if not any(_is_overlapped(donor_start, donor_end, left, right) for left, right in blocks):
            continue
        # Non-primary records cannot override the disposition of a valid pair.
        if rec.flag & (0x4 | 0x100 | 0x800):
            continue
        if not rec.qname or rec.qname == "*":
            raise ValueError("a donor-overlapping primary alignment needs a query name to identify its fragment")
        key = fragment_key(rec)
        if not _primary_sam(rec, min_mapq):
            res.overlapped.add(key)
            res.discarded.add(key)
            continue
        res.overlapped.add(key)
        res.flag[key] = rec.flag
        duplicate = decode_duplicate(rec.flag)
        if (not rec.flag & 0x1 or rec.flag & 0x8 or rec.rnext == "*" or rec.pnext < 1
                or bool(rec.flag & 0x40) == bool(rec.flag & 0x80)
                or (rec.rnext in ("=", rec.rname) and rec.tlen == 0)):
            res.discarded.add(key)
            continue
        near = rec.rnext in ("=", rec.rname) and 0 < abs(rec.tlen) < 1000
        if near:
            frag_start = min(rec.pos, rec.pnext)
            frag_end = frag_start + abs(rec.tlen) - 1
            if rec.mate_cigar and rec.mate_cigar != "*":
                # planning/328: in a DOVETAILED short fragment the reverse mate starts left of the forward mate,
                # and bwa's TLEN (measured from the forward mate) does not span the reverse read's alignment, so
                # the TLEN test below discarded a genomic fragment that runs past the donor into the chromosome.
                # With the mate's CIGAR the fragment's aligned extent is known exactly.
                frag_end = max(read_end, _cigar_blocks(rec.pnext, rec.mate_cigar)[1])
            # TLEN must actually contain the represented read and mate start.
            if frag_end < max(read_end, rec.pnext):
                res.discarded.add(key)
                continue
            template = _is_contained(frag_start, frag_end, donor_start - flexible_bases, donor_end + flexible_bases)
            chimeric = ((leading_soft and abs(rec.pos - donor_start) <= flexible_bases)
                        or (trailing_soft and abs(read_end - donor_end) <= flexible_bases))
        else:
            mate_chr = rec.rname if rec.rnext == "=" else rec.rnext
            cassette_mate = (mate_chr in plasmid_contigs
                             or any(start <= rec.pnext <= end for start, end in blacklist.get(mate_chr, ())))
            template = cassette_mate and _is_contained(rec.pos, read_end, donor_start - flexible_bases, donor_end + flexible_bases)
            chimeric = cassette_mate and not template
        if template and keep_contained and near and not chimeric:
            # Plasmids purged: a contained same-chromosome fragment can only be genomic.
            res.included.add(key)
            res.contained_kept.add(key)
            if not duplicate:
                res.unique_contained_kept.add(key)
        elif template:
            res.discarded.add(key)
            if not duplicate:
                res.unique_plasmid.add(key)
        else:
            res.included.add(key)
            if chimeric:
                res.redi_translocation.add(key)
                if not duplicate:
                    res.unique_chimeric.add(key)
            elif not duplicate:
                res.unique_genome.add(key)
                if any(left <= midpoint <= right for left, right in blocks):
                    coverage[key] = coverage.get(key, 0) + 1
    res.included -= res.discarded
    res.redi_translocation -= res.discarded
    res.unique_chimeric -= res.discarded
    res.contained_kept -= res.discarded
    res.unique_contained_kept -= res.discarded
    res.unique_genome -= res.discarded | res.unique_chimeric | res.contained_kept
    res.variant_cov = sum(coverage.get(key, 0) for key in res.unique_genome)
    return res


# --------------------------------------------------------------------------- #
# Plasmid load: when may donor-contained fragments be kept? (planning/328)     #
# --------------------------------------------------------------------------- #
#
# Kevin R. Roy, 2026-10-30: "count contained read pairs when the backbone is empty. This is safe when we know
# that we have removed plasmids, so it needs to be a per clone/per dataset option for trace, that when it is
# high confidence that the plasmids have been purged, we don't need to restrict to segments that span the
# donor-genome junction." A guide-donor plasmid present at k copies per genome leaves reads across its
# backbone at k times nuclear depth and pairs a donor read with a backbone mate; both are measured here.
# Measured on the NNS cohort (planning/328s, 328t): 24 backbone reads over 523 clones, no donor-backbone pair.

_NUCLEAR_RE = re.compile(r"chr[IVX]+")
Interval = Tuple[str, int, int]


def poisson_upper(n: int, tail: float = 0.025) -> float:
    """Exact upper confidence limit for a Poisson mean given ``n`` observed events.

    The lambda at which P(X <= n) equals ``tail`` (0.025 gives the one-sided 97.5% limit, 3.689 for n = 0),
    solved by bisection so the module needs no scipy.
    """
    import math

    if type(n) is not int or n < 0 or not 0 < tail < 1:
        raise ValueError("n must be a non-negative int and tail in (0, 1)")

    def cdf(lam: float) -> float:
        return sum(math.exp(k * math.log(lam) - lam - math.lgamma(k + 1)) for k in range(n + 1)) if lam > 0 else 1.0

    lo, hi = 0.0, n + 20.0 * math.sqrt(n + 1) + 20.0
    for _ in range(200):
        mid = (lo + hi) / 2
        if cdf(mid) > tail:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


@dataclass(frozen=True)
class PlasmidLoad:
    """Evidence for or against a plasmid in one sample.

    ``plasmid_reads`` are primary, non-QC-fail reads with MAPQ >= ``min_mapq`` on the declared plasmid contigs,
    outside the masked intervals (duplicates included, as in the index statistics the nuclear count comes
    from). ``donor_backbone_pairs`` counts fragments with a read inside a donor window +/- slack and its mate on
    a declared plasmid contig. ``copies`` is plasmid read density over nuclear read density; ``copies_upper95``
    its exact Poisson 97.5% upper limit, infinite when nuclear depth is zero.
    """

    plasmid_reads: int
    unmasked_length: int
    nuclear_reads: int
    nuclear_length: int
    donor_backbone_pairs: int
    expected_reads_per_copy: float
    copies: float
    copies_upper95: float


def plasmid_load(
    bam_path: str,
    plasmid_contigs: Sequence[str],
    *,
    mask: Iterable[Interval] = (),
    donor_windows: Iterable[Interval] = (),
    nuclear_contigs: Optional[Sequence[str]] = None,
    min_mapq: int = 1,
    flexible_bases: int = FLEXIBLE_BASES,
) -> PlasmidLoad:
    """Measure plasmid load in one indexed BAM (see :class:`PlasmidLoad`).

    ``plasmid_contigs`` should be the plasmids that CARRY the donor; a plasmid without donor sequence cannot make
    a donor-contained fragment and only dilutes the test. ``mask`` excludes 1-based inclusive intervals that
    collect reads without any plasmid present, such as a bc1 cassette that amplicon product maps to.
    ``nuclear_contigs`` defaults to every ``chr<roman>`` contig in the header.
    """
    import pysam

    plasmid_contigs = list(plasmid_contigs)
    if not plasmid_contigs:
        raise ValueError("plasmid_load needs at least one plasmid contig")
    if type(min_mapq) is not int or not 0 <= min_mapq <= 254:
        raise ValueError("min_mapq must be an int from 0 to 254")
    mask_by: Dict[str, List[Tuple[int, int]]] = {}
    for contig, start, end in mask:
        if start < 1 or end < start:
            raise ValueError(f"invalid mask interval {contig}:{start}-{end}")
        mask_by.setdefault(contig, []).append((start, end))
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        lengths = dict(zip(bam.references, bam.lengths))
        missing = [c for c in plasmid_contigs if c not in lengths]
        if missing:
            raise ValueError(f"plasmid contig(s) absent from the BAM header: {missing}")
        for contig, spans in mask_by.items():
            if contig not in lengths or any(end > lengths[contig] for _, end in spans):
                raise ValueError(f"mask interval absent or outside the BAM reference: {contig}")
        nuclear = list(nuclear_contigs) if nuclear_contigs is not None else [c for c in lengths if _NUCLEAR_RE.fullmatch(c)]
        if not nuclear or any(c not in lengths for c in nuclear):
            raise ValueError("nuclear contigs are missing from the BAM header; pass nuclear_contigs explicitly")
        stats = {s.contig: s.mapped for s in bam.get_index_statistics()}
        nuclear_reads = sum(stats.get(c, 0) for c in nuclear)
        nuclear_length = sum(lengths[c] for c in nuclear)

        plasmid_reads = 0
        masked_length = 0
        for contig in plasmid_contigs:
            spans = mask_by.get(contig, [])
            covered = set()
            for start, end in spans:
                covered.update(range(start, end + 1))
            masked_length += len(covered)
            for read in bam.fetch(contig):
                if read.flag & (0x4 | 0x100 | 0x200 | 0x800) or read.mapping_quality < min_mapq:
                    continue
                left, right = read.reference_start + 1, read.reference_end
                if any(left <= end and start <= right for start, end in spans):
                    continue
                plasmid_reads += 1

        pairs = set()
        declared = set(plasmid_contigs)
        for chrom, start, end in donor_windows:
            if chrom not in lengths:
                raise ValueError(f"donor window contig absent from the BAM header: {chrom}")
            for read in bam.fetch(chrom, max(0, start - 1 - flexible_bases), end + flexible_bases):
                if read.flag & (0x4 | 0x100 | 0x200 | 0x800) or read.reference_end is None:
                    continue
                if not read.is_paired or read.mate_is_unmapped or read.next_reference_name not in declared:
                    continue
                if read.reference_start + 1 >= start - flexible_bases and read.reference_end <= end + flexible_bases:
                    pairs.add((read.get_tag("RG") if read.has_tag("RG") else "", read.query_name))

    unmasked = sum(lengths[c] for c in plasmid_contigs) - masked_length
    per_bp = nuclear_reads / nuclear_length if nuclear_length else 0.0
    expected = per_bp * unmasked
    copies = plasmid_reads / expected if expected else float("inf")
    upper = poisson_upper(plasmid_reads) / expected if expected else float("inf")
    return PlasmidLoad(plasmid_reads, unmasked, nuclear_reads, nuclear_length, len(pairs), expected, copies, upper)


def plasmids_purged(load: PlasmidLoad, max_copies: float = 0.05) -> bool:
    """High-confidence purge: no donor-backbone pair and the copy number's upper limit at or below ``max_copies``.

    The default 0.05 copies per genome bounds the share of donor-contained fragments that could be template near
    5%, the same order as the foreign-genome share measured at the bc1 locus (planning/328u). A sample too shallow
    to bound its copy number fails the test and keeps the historical discard.
    """
    if not max_copies > 0:
        raise ValueError("max_copies must be positive")
    return load.donor_backbone_pairs == 0 and load.copies_upper95 <= max_copies


def _iter_sam(bam_path: str, samtools: str = "samtools") -> Iterator[str]:
    """Yield SAM lines (incl. ``@`` header) from ``samtools view -h <bam>``.

    Raises :class:`subprocess.CalledProcessError` if ``samtools`` exits non-zero
    (e.g. a truncated or unreadable BAM), so callers never silently process a
    partial SAM stream and emit wrong counts. The exit code is only checked on
    normal exhaustion of the generator; a consumer that stops early triggers
    ``GeneratorExit`` and skips the check (no spurious error on a broken pipe).
    """
    proc = subprocess.Popen(
        [samtools, "view", "--no-PG", "-h", bam_path],
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


_REFERENCE_SIGNATURES: Dict[tuple, Dict[str, tuple[int, str]]] = {}


def _reference_signature(path: str) -> Dict[str, tuple[int, str]]:
    """Sequence identity by named contig, cached across samples of one cohort."""
    import pysam
    stat, fai = os.stat(path), os.stat(path + ".fai")
    key = (os.path.realpath(path), stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns,
           fai.st_size, fai.st_mtime_ns, fai.st_ctime_ns)
    if key not in _REFERENCE_SIGNATURES:
        result = {}
        with pysam.FastaFile(path) as fasta:
            for contig, length in zip(fasta.references, fasta.lengths):
                digest = hashlib.sha256()
                for start in range(0, length, 1048576):
                    digest.update(fasta.fetch(contig, start, min(length, start + 1048576)).upper().encode("ascii"))
                result[contig] = (length, digest.hexdigest())
        _REFERENCE_SIGNATURES[key] = result
    return _REFERENCE_SIGNATURES[key]


def bam_index_suffix(reference_lengths: Iterable[int]) -> str:
    """BAI for ordinary contigs, CSI when a reference reaches BAI's 2^29 limit."""
    return ".csi" if any(length >= (1 << 29) for length in reference_lengths) else ".bai"


def _remove_bam_indexes(path: str) -> None:
    """Remove every htslib-discoverable old index when replacing its BAM."""
    for prefix in (path, os.path.splitext(path)[0]):
        for suffix in (".bai", ".csi"):
            try:
                os.unlink(prefix + suffix)
            except FileNotFoundError:
                pass


@contextmanager
def output_lock(path: str):
    """Serialize writers of one output bundle using a persistent advisory lock.

    Lock files are intentionally retained: unlinking a released lock can let an
    existing waiter and a new caller lock different inodes. A process exit releases
    the OS lock automatically; file existence does not indicate an active writer.
    All writers must honor this protocol. Supported execution hosts are POSIX.
    """
    import fcntl
    parent = os.path.realpath(os.path.dirname(os.path.abspath(path)))
    os.makedirs(parent, exist_ok=True)
    lock_path = os.path.join(parent, os.path.basename(path) + ".lock")
    with open(lock_path, "a") as handle:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


def clean_donor_reads(
    in_bam: str, out_bam: str, chrom: str, donor_start: int, donor_end: int,
    blacklist_path: str, flexible_bases: int = FLEXIBLE_BASES,
    dup_mode: str = "standard", samtools: str = "samtools", min_mapq: int = 20,
    reference_path: Optional[str] = None, blacklist_reference: Optional[str] = None,
    *, keep_contained: bool = False, plasmid_contigs: Iterable[str] = (),
) -> CleanResult:
    """Publish a filtered, indexed full BAM and fragment annotation under one lock.

    Standard mode requires identical declared sample/blacklist reference sequences.
    Filtering preserves coordinate order and original headers. Conversion, indexing
    and validation complete before publication; concurrent callers for the same
    output wait for the current generation. Explicit ``dup_mode="parity"`` retains
    historical counting semantics. The annotation describes retained genome,
    discarded donor-template and candidate chimeric fragments at one donor window.
    """
    with output_lock(out_bam):
        return _clean_donor_reads_unlocked(
            in_bam, out_bam, chrom, donor_start, donor_end, blacklist_path,
            flexible_bases, dup_mode, samtools, min_mapq, reference_path, blacklist_reference,
            keep_contained=keep_contained, plasmid_contigs=plasmid_contigs,
        )


def _clean_donor_reads_unlocked(
    in_bam: str,
    out_bam: str,
    chrom: str,
    donor_start: int,
    donor_end: int,
    blacklist_path: str,
    flexible_bases: int = FLEXIBLE_BASES,
    dup_mode: str = "standard",
    samtools: str = "samtools",
    min_mapq: int = 20,
    reference_path: Optional[str] = None,
    blacklist_reference: Optional[str] = None,
    *,
    keep_contained: bool = False,
    plasmid_contigs: Iterable[str] = (),
) -> CleanResult:
    """Filter donor/plasmid/unresolved reads from a target site (Half B driver).

    Two passes over the input BAM (exactly like the Perl):

    1. classify every donor-overlapping read pair and write
       ``{out_bam}.target_filtering.anno``;
    2. write ``out_bam`` keeping all non-target reads plus the donor-overlapping
       reads that were *included* and not *discarded* (**discard beats include**).

    Standard mode additionally requires ``reference_path`` and
    ``blacklist_reference`` to have identical named-contig sequences. An explicitly
    declared renamed FASTA is accepted; coordinates are never lifted by name.
    Mode affects both evidence eligibility and cleaned-BAM membership. The
    BAM, index and annotation are published only after conversion and indexing
    succeed. Input coordinate order is preserved; an unindexable unsorted input
    fails rather than becoming a success-looking cleaned artifact. Parity is for historical
    replication and is not the production filtering policy.

    Returns the :class:`CleanResult` from pass 1.
    """
    import pysam

    if shutil.which(samtools) is None:
        raise FileNotFoundError(f"samtools executable not found: {samtools}")
    if os.path.realpath(in_bam) == os.path.realpath(out_bam):
        raise ValueError("cleaned BAM must differ from input BAM")
    blacklist = load_blacklist(blacklist_path, mode=dup_mode)
    with pysam.AlignmentFile(in_bam, "rb") as bam:
        lengths = dict(zip(bam.references, bam.lengths))
    if chrom not in lengths or not 1 <= donor_start <= donor_end <= lengths[chrom]:
        raise ValueError(f"donor interval is absent or outside BAM reference: {chrom}:{donor_start}-{donor_end}")
    if dup_mode != "parity":
        if reference_path is None or blacklist_reference is None:
            raise ValueError("standard donor cleaning requires reference_path and blacklist_reference provenance")
        reference_signature = _reference_signature(reference_path)
        if {ref: signature[0] for ref, signature in reference_signature.items()} != lengths:
            raise ValueError("BAM and supplied reference contig dictionaries disagree")
        if reference_signature != _reference_signature(blacklist_reference):
            raise ValueError("blacklist reference sequence identity differs from sample reference; coordinate lifting is required")
        if not blacklist:
            raise ValueError("blacklist contains no reference intervals")
        for ref, intervals in blacklist.items():
            if ref not in lengths or any(end > lengths[ref] for _, end in intervals):
                raise ValueError(f"blacklist interval absent or outside BAM reference: {ref}")

    def records():
        for line in _iter_sam(in_bam, samtools):
            if line.startswith("@"):
                continue
            rec = parse_sam_line(line)
            mate_chr = rec.rname if rec.rnext == "=" else rec.rnext
            if dup_mode != "parity" and (mate_chr not in lengths or rec.pnext > lengths[mate_chr]):
                rec.flag |= 0x8
            yield rec

    result = classify_donor_reads(records(), chrom, donor_start, donor_end, blacklist, flexible_bases, dup_mode, min_mapq,
                                  keep_contained=keep_contained, plasmid_contigs=plasmid_contigs)
    # Pass 2 must look records up in the namespace pass 1 recorded them under.
    # The bare-name prefilter only skips building a key for untouched records.
    if dup_mode == "parity":
        identity = lambda rec: rec.qname  # noqa: E731 - historical QNAME-only namespace
        touched_names = result.overlapped
    else:
        identity = fragment_key
        touched_names = {key.qname for key in result.overlapped}
    parent = os.path.dirname(os.path.abspath(out_bam))
    os.makedirs(parent, exist_ok=True)
    # Publish only after both samtools passes, conversion and BAM validation have
    # succeeded. A failed subprocess cannot leave a success-looking .anno file.
    with tempfile.TemporaryDirectory(prefix=".donor_clean_", dir=parent) as tmp:
        tmp_sam = os.path.join(tmp, "clean.sam")
        tmp_bam = os.path.join(tmp, "clean.bam")
        tmp_anno = os.path.join(tmp, "clean.anno")
        index_suffix = bam_index_suffix(lengths.values())
        tmp_index = tmp_bam + index_suffix
        with open(tmp_sam, "w") as out:
            for line in _iter_sam(in_bam, samtools):
                if line.startswith("@"):
                    out.write(line)
                    continue
                rec = parse_sam_line(line)
                if dup_mode != "parity" and rec.rname == chrom and rec.pos > 0:
                    blocks, _, _, _ = _cigar_blocks(rec.pos, rec.cigar)
                    at_target = any(_is_overlapped(donor_start, donor_end, left, right) for left, right in blocks)
                    if at_target and not _primary_sam(rec, min_mapq):
                        continue
                if rec.qname not in touched_names:
                    out.write(line)
                    continue
                key = identity(rec)
                if key not in result.overlapped or (key in result.included and key not in result.discarded):
                    out.write(line)
        with open(tmp_bam, "wb") as out_fh:
            # Temporary paths must not become nondeterministic @PG header data.
            # The stage cache/provenance records filtering inputs and source.
            subprocess.run([samtools, "view", "--no-PG", "-bS", tmp_sam], stdout=out_fh, check=True)
        with pysam.AlignmentFile(tmp_bam, "rb") as checked:
            for _ in checked:
                pass
        pysam.index(*(("-c", tmp_bam) if index_suffix == ".csi" else (tmp_bam,)))
        with pysam.AlignmentFile(tmp_bam, "rb") as checked:
            checked.check_index()
        write_anno(result, tmp_anno)
        _remove_bam_indexes(out_bam)
        os.replace(tmp_bam, out_bam)
        os.replace(tmp_index, out_bam + index_suffix)
        os.replace(tmp_anno, out_bam + ".target_filtering.anno")
    return result


# --------------------------------------------------------------------------- #
# Half A: plasmid-integration detection                                        #
# --------------------------------------------------------------------------- #


def _integrated_contigs(references: Sequence[str]) -> List[str]:
    return [r for r in references if r.startswith("integrated_plasmid_")]


def detect_plasmid_integration(
    bam_path: str,
    integrated_contigs: Optional[Sequence[str]] = None,
    junction_distance: Optional[Union[int, Dict[str, int]]] = None,
    window: int = 200,
    min_mapq: int = 20,
    min_anchor: int = 1,
) -> Dict[str, int]:
    """Count primary reads spanning the genome↔integration junctions of a sample.

    Requires aligned reference bases on both sides of a known junction. D/N spans
    and reads starting exactly at a junction do not qualify. Primary, QC-passing
    records with MAPQ >= ``min_mapq`` count; duplicates remain raw read records.
    MAPQ 20 and one aligned base per side are policy defaults, not calibrated
    integration sensitivity. A per-contig distance must be supplied, derived as
    ``upstream_arm_length + donor_length``. An explicit scalar is allowed for
    known uniform layouts; there is no implicit 1129 assumption. Missing target
    contigs, missing distances and zero eligible contigs raise, never return a
    biological-negative zero for an unmeasured target.

    Returns ``{"upstream_count", "downstream_count"}`` summed over all contigs for
    this (single-sample) BAM.
    """
    import pysam

    if (any(type(value) is not int for value in (window, min_anchor, min_mapq))
            or window <= 0 or min_anchor < 1 or not 0 <= min_mapq <= 254):
        raise ValueError("invalid junction window, minimum anchor or MAPQ threshold")
    bam = pysam.AlignmentFile(bam_path, "rb")
    try:
        refs = bam.references
        contigs = list(integrated_contigs) if integrated_contigs is not None else _integrated_contigs(refs)
        if not contigs:
            raise ValueError("no eligible integrated-plasmid contigs; integration is unmeasured")
        if len(set(contigs)) != len(contigs):
            raise ValueError("duplicate integrated-plasmid contigs")
        distances: Dict[str, int] = {}
        for ref in contigs:
            if ref not in refs:
                raise ValueError(f"integrated contig absent from BAM reference: {ref}")
            jd = junction_distance.get(ref) if isinstance(junction_distance, dict) else junction_distance
            length = bam.get_reference_length(ref)
            if isinstance(jd, bool) or not isinstance(jd, int) or not min_anchor <= jd < length - jd:
                raise ValueError(f"missing or invalid junction distance for {ref}: {jd!r} (length {length})")
            if jd + min_anchor > length - jd:
                raise ValueError(f"junction anchors overlap on {ref}")
            distances[ref] = jd
        up = 0
        down = 0
        for ref in contigs:
            jd = distances[ref]
            up_junction = jd
            down_junction = bam.get_reference_length(ref) - jd
            for junction, side in ((up_junction, "up"), (down_junction, "down")):
                start = max(0, junction - window)
                end = min(bam.get_reference_length(ref), junction + window)
                for read in bam.fetch(ref, start, end):
                    if not _eligible_sv_read(read, min_mapq):
                        continue
                    # Merge abutting M/= /X blocks, including those separated by
                    # an insertion; a D/N leaves a reference gap and cannot merge.
                    blocks = _merged_aligned_blocks(read)
                    if any(left <= junction - min_anchor and right >= junction + min_anchor for left, right in blocks):
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
    upstream_count: Optional[int]
    downstream_count: Optional[int]
    total_count: Optional[int]
    half_a_status: str


def summarize_plasmid_integration(
    sample_bams: Dict[str, str],
    integrated_contigs: Optional[Sequence[str]] = None,
    junction_distance: Optional[Union[int, Dict[str, int]]] = None,
    window: int = 200,
) -> List[PlasmidIntegrationRow]:
    """Run :func:`detect_plasmid_integration` over ``{sample_name: bam_path}``.

    Returns one :class:`PlasmidIntegrationRow` per sample, sorted by ``total_count``
    descending — the row schema of Kevin's ``plasmid_integration_summary.tsv``.
    """
    rows: List[PlasmidIntegrationRow] = []
    for sample, bam_path in sample_bams.items():
        import pysam
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            contigs = list(integrated_contigs) if integrated_contigs is not None else _integrated_contigs(bam.references)
        if not contigs:
            rows.append({"sample_name": sample, "upstream_count": None, "downstream_count": None,
                         "total_count": None, "half_a_status": "not_run_no_integrated_contigs"})
            continue
        counts = detect_plasmid_integration(bam_path, contigs, junction_distance, window)
        total = counts["upstream_count"] + counts["downstream_count"]
        rows.append(
            {
                "sample_name": sample,
                "upstream_count": counts["upstream_count"],
                "downstream_count": counts["downstream_count"],
                "total_count": total,
                "half_a_status": "run",
            }
        )
    rows.sort(key=lambda r: -1 if r["total_count"] is None else r["total_count"], reverse=True)
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
    if type(min_gap) is not int or min_gap < 0:
        raise ValueError("min_gap must be nonnegative")
    if min(r1_start, r2_start) < 0 or r1_end < r1_start or r2_end < r2_start:
        return False
    region_length = max(r1_end, r2_end) - min(r1_start, r2_start) + 1
    r1_length = r1_end - r1_start + 1
    r2_length = r2_end - r2_start + 1
    if region_length > r2_length + min_gap and (r2_start > r1_start and r2_end < r1_end):
        return True
    if region_length > r1_length + min_gap and (r1_start > r2_start and r1_end < r2_end):
        return True
    return False


def _eligible_sv_read(read, min_mapq: int = 20) -> bool:
    return (not (read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_qcfail)
            and min_mapq <= read.mapping_quality < 255 and read.reference_id >= 0)


def _merged_aligned_blocks(read) -> list[tuple[int, int]]:
    blocks: list[tuple[int, int]] = []
    for left, right in read.get_blocks():
        if blocks and blocks[-1][1] == left:
            blocks[-1] = (blocks[-1][0], right)
        else:
            blocks.append((left, right))
    return blocks


def mate_pair_large_deletion(read_1, read_2, min_gap: int = 1000, min_mapq: int = 20) -> bool:
    """pysam wrapper for :func:`one_read_inside_large_deletion`.

    Uses each read's first/last aligned reference position
    (``get_reference_positions()[0]`` / ``[-1]``), matching the original.
    """
    if any(type(v) is not int for v in (min_gap, min_mapq)) or min_gap < 0 or not 0 <= min_mapq <= 254:
        raise ValueError("invalid deletion gap or MAPQ threshold")
    if not all(_eligible_sv_read(r, min_mapq) and r.is_paired and not r.mate_is_unmapped for r in (read_1, read_2)):
        return False
    if any(op == 3 for read in (read_1, read_2) for op, _ in (read.cigartuples or ())):
        return False
    if (not read_1.query_name or read_1.query_name == "*" or read_1.query_name != read_2.query_name
            or read_1.reference_id != read_2.reference_id
            or read_1.reference_name != read_2.reference_name
            or not ((read_1.is_read1 and read_2.is_read2) or (read_1.is_read2 and read_2.is_read1))
            or read_1.is_read1 == read_1.is_read2 or read_2.is_read1 == read_2.is_read2
            or read_1.next_reference_id != read_2.reference_id
            or read_2.next_reference_id != read_1.reference_id
            or read_1.next_reference_start != read_2.reference_start
            or read_2.next_reference_start != read_1.reference_start):
        return False
    p1 = read_1.get_blocks()
    p2 = read_2.get_blocks()
    if not p1 or not p2:
        return False
    return one_read_inside_large_deletion(p1[0][0], p1[-1][1] - 1, p2[0][0], p2[-1][1] - 1, min_gap)


def within_read_large_deletions(read, cutoff: int = 50, min_mapq: int = 20):
    """Within-read deletions strictly larger than ``cutoff`` bp.

    Reuses :func:`trace_crispr.core.cigar.get_deletions_from_cigar`. The original
    classifies a deletion as "large" when ``len(deleted_bases) > cutoff`` (strict),
    i.e. size ``>= cutoff + 1``.
    """
    if any(type(v) is not int for v in (cutoff, min_mapq)) or cutoff < 0 or not 0 <= min_mapq <= 254:
        raise ValueError("invalid deletion cutoff or MAPQ threshold")
    if not _eligible_sv_read(read, min_mapq):
        return []
    blocks = read.get_blocks()
    return [deletion for deletion in get_deletions_from_cigar(read, min_size=cutoff + 1)
            if any(left <= deletion.ref_start - 1 < right for left, right in blocks)
            and any(left <= deletion.ref_end < right for left, right in blocks)]


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
    if not plasmid_sequence or type(backbone_length) is not int or not 0 <= backbone_length <= len(plasmid_sequence):
        raise ValueError("backbone length must lie within a nonempty plasmid sequence")
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
    if donor_strand not in ("+", "-"):
        raise ValueError("donor_strand must be '+' or '-'")
    if type(bp_to_extend) is not int or bp_to_extend < 1 or not donor or not integrated_segment or len(extended_donor) < 2 * bp_to_extend + len(donor):
        raise ValueError("integrated contig requires two complete positive-length arms, donor and inserted segment")
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
    if type(bp_to_extend) is not int or bp_to_extend < 1 or not donor or not re.fullmatch(r"[ACGTRYSWKMBDHVNacgtryswkmbdhvn]+", donor):
        raise ValueError("junction derivation requires positive arm length and a nonempty DNA donor")
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
    p.add_argument("--min-mapq", type=int, default=20)
    p.add_argument("--reference", dest="reference_path", help="sample alignment reference FASTA (required in standard mode)")
    p.add_argument("--blacklist-reference", help="source reference FASTA for blacklist coordinates (required in standard mode)")
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
        min_mapq=args.min_mapq,
        reference_path=args.reference_path,
        blacklist_reference=args.blacklist_reference,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
