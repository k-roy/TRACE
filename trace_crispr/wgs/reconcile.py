"""
trace_crispr.wgs.reconcile
==========================

WGS Stage 3 — **reconciliation**. Combine the per-cell bc1 lineage set
(:mod:`~trace_crispr.wgs.bc1_consensus`) with the per-designed-site edit VAF to
classify each WGS sample as ``NO_EDIT`` / ``SINGLE_EDIT`` / ``MULTI_EDIT``
(a genuine cell carrying >=2 integrated donors) / ``CONTAMINATION`` (a mixture of
lineages), with a transparent confidence label.

Biology (haploid yeast, single-colony picks) — validated 2026-06-04:

* A clone's edit at a designed locus is **bimodal**: VAF~0 (WT) or VAF~1 (edited).
  Across 450 pure rearrayed wells only 0.4% were intermediate, so an **intermediate
  VAF at a designed site => a mixed founder => contamination** (mosaic / CNV /
  incomplete editing do not create an intermediate-VAF population).
* The two WGS bc1 stages == the two extraction tracks. The **plasmid track**
  (pre-integration) sees plasmid-borne bc1, so a cell that took up two distinct
  guide-donor plasmids shows >=2 bc1 there — *genuine*, not contamination. The
  **genome track** (post-integration, yT177 integrant) has a single haploid locus
  that only one plasmid wins, so **>=2 bc1 on the genome track => a mixture**.
* A genuine ``MULTI_EDIT`` is therefore two plasmids / two oligos -> **two genomic
  edits, both at VAF~1** (both donors repaired the one haploid genome), NOT one
  oligo encoding two edits.

Design consequence (the central architectural choice): the **edit count is read
from genomic VAF at designed bases** (robust, bridge-free), while bc1 bridging is
used only for *lineage attribution* and the *genome-track >=2 = mixture* signal.
This keeps the genuine-double determination from being held hostage to bridging
both barcodes (the failure mode of the first PTC verification pass).

The classification core (:func:`reconcile_cell`, :func:`classify_site`) is pure and
takes already-measured site VAFs, so it is fully unit-testable without a BAM; the
I/O layer (:func:`window_edit_vaf`, :func:`build_edit_sites`, :func:`run`) wires it
to a BAM + a per-library ``oligo -> donor window`` map.
"""

from __future__ import annotations

import argparse
import csv
import sys
from dataclasses import dataclass, field
from difflib import SequenceMatcher
from typing import Dict, Iterable, List, Optional, Tuple

try:  # pysam is only needed for the BAM I/O layer, not the pure core/tests
    import pysam
except ImportError:  # pragma: no cover
    pysam = None  # type: ignore

# --------------------------------------------------------------------------- #
# Labels                                                                       #
# --------------------------------------------------------------------------- #

# Per-site status
SITE_FULL = "full"  # VAF >= tau_full, covered            -> an edit is present
SITE_PARTIAL = "partial"  # tau_wt < VAF < tau_full, covered -> mixture signature
SITE_WT = "wt"  # VAF <= tau_wt, covered                  -> not edited here
SITE_NOCOV = "nocov"  # depth < min_cov                    -> cannot confirm

# Per-cell label
NO_EDIT = "NO_EDIT"
SINGLE_EDIT = "SINGLE_EDIT"
MULTI_EDIT = "MULTI_EDIT"
CONTAMINATION = "CONTAMINATION"
LOW_CONF = "LOW_CONF"  # nothing covered well enough to classify


# --------------------------------------------------------------------------- #
# Data                                                                         #
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class EditSite:
    """A designed edit window to score, with optional lineage attribution.

    ``start``/``end`` are 1-based inclusive genomic coordinates of the donor window.
    ``oligo`` records which bridged bc1->oligo predicted this site (None when the
    site comes from a library-wide designed-site scan rather than this cell's bc1).
    ``wt_seq``/``design_donor`` are optional but strongly preferred: when present,
    VAF is scored only at designed bases rather than every non-reference base in
    the donor window.
    """

    chrom: str
    start: int
    end: int
    oligo: Optional[str] = None
    wt_seq: Optional[str] = None
    design_donor: Optional[str] = None

    def key(self) -> Tuple[str, int, int, Optional[str], Optional[str]]:
        """Stable identity for de-duplicating library rows without collapsing
        distinct designs that happen to share a donor window."""
        return (self.chrom, self.start, self.end, self.wt_seq, self.design_donor)


@dataclass
class SiteCall:
    """One scored designed site within a cell."""

    site: EditSite
    vaf: float
    cov: int
    status: str  # SITE_FULL / SITE_PARTIAL / SITE_WT / SITE_NOCOV


@dataclass
class CellReconciliation:
    """The Stage-3 classification of one WGS sample."""

    sample: str
    label: str
    n_genome_bc: int
    n_plasmid_bc: int
    n_sites: int
    n_covered: int
    n_full: int
    n_partial: int
    n_wt: int
    n_nocov: int
    sites: List[SiteCall] = field(default_factory=list)
    confidence: str = "low"  # high / low
    flags: Tuple[str, ...] = ()


# --------------------------------------------------------------------------- #
# Pure classification core (no BAM — fully unit-testable)                      #
# --------------------------------------------------------------------------- #


def classify_site(
    vaf: float,
    cov: int,
    *,
    tau_full: float = 0.85,
    tau_wt: float = 0.15,
    min_cov: int = 4,
    min_partial_reads: int = 3,
) -> str:
    """Bucket one designed site by its edit VAF and coverage.

    ``min_partial_reads`` is an **absolute non-reference read-count floor** on the
    PARTIAL (mixture) call, kept deliberately separate from the VAF thresholds. At
    ``min_cov``=4 a single spurious non-ref read is VAF 0.25 — in the PARTIAL band —
    so scanning ~485 designed sites/cell would inflate false CONTAMINATION. A genuine
    minor-lineage secondary is ~30% of a full-depth well (~18 of 60 reads); the
    absolute count (18 vs 1) separates real contamination from noise on a *single*
    site, without raising ``tau_wt`` (which would discard the genuine low-VAF
    secondary edits we are trying to recover). Below the floor, an intermediate VAF
    is treated as WT rather than PARTIAL.
    """
    if cov < min_cov:
        return SITE_NOCOV
    if vaf >= tau_full:
        return SITE_FULL
    if vaf <= tau_wt:
        return SITE_WT
    if round(vaf * cov) < min_partial_reads:
        return SITE_WT
    return SITE_PARTIAL


def reconcile_cell(
    sample: str,
    sites: List[SiteCall],
    *,
    n_genome_bc: int = 0,
    n_plasmid_bc: int = 0,
    bc_confidence: str = "high",
) -> CellReconciliation:
    """Classify a cell from its already-scored designed sites + bc1 multiplicity.

    Decision order (a mixture signal dominates a clean multi-edit):
      1. any PARTIAL site, or >=2 GENOME-track bc1   -> CONTAMINATION
      2. >=2 FULL sites                              -> MULTI_EDIT
      3. exactly 1 FULL site                         -> SINGLE_EDIT
      4. >=1 covered (WT) site, 0 full/partial       -> NO_EDIT
      5. nothing covered                             -> LOW_CONF

    Rationale: a haploid clone is 0/100 at each locus, so an intermediate (PARTIAL)
    VAF can only come from a mixed founder; and only one plasmid integrates the
    genome locus, so >=2 genome-track bc1 is itself a mixture — both override the
    "two clean edits" reading.
    """
    n_sites = len(sites)
    n_full = sum(1 for s in sites if s.status == SITE_FULL)
    n_partial = sum(1 for s in sites if s.status == SITE_PARTIAL)
    n_wt = sum(1 for s in sites if s.status == SITE_WT)
    n_nocov = sum(1 for s in sites if s.status == SITE_NOCOV)
    n_covered = n_full + n_partial + n_wt

    flags: List[str] = []
    if n_partial:
        flags.append("PARTIAL_VAF")
    if n_genome_bc >= 2:
        flags.append("MULTI_GENOME_BC")

    if n_partial >= 1 or n_genome_bc >= 2:
        label = CONTAMINATION
    elif n_full >= 2:
        label = MULTI_EDIT
    elif n_full == 1:
        label = SINGLE_EDIT
    elif n_covered >= 1:
        label = NO_EDIT
    else:
        label = LOW_CONF

    # Confidence: trust the call only when the bc1 call is high-confidence and no
    # designed site needed to rule out an edit is uncovered. A NO_EDIT resting on
    # uncovered sites is not trustworthy (an edit there could simply be unsequenced).
    confidence = "high" if bc_confidence == "high" else "low"
    if label == LOW_CONF:
        confidence = "low"
    elif label == NO_EDIT and n_nocov:
        confidence = "low"
        flags.append("UNCOVERED_SITES")

    return CellReconciliation(
        sample=sample, label=label, n_genome_bc=n_genome_bc, n_plasmid_bc=n_plasmid_bc,
        n_sites=n_sites, n_covered=n_covered, n_full=n_full, n_partial=n_partial,
        n_wt=n_wt, n_nocov=n_nocov, sites=sites,
        confidence=confidence, flags=tuple(flags),
    )


# --------------------------------------------------------------------------- #
# BAM I/O layer                                                                #
# --------------------------------------------------------------------------- #

_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _rev_comp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1].upper()


def _designed_offsets(wt_seq: str, design_donor: str) -> List[int]:
    """Return 0-based reference offsets affected by the designed donor edit.

    Substitutions and deletions mark the changed/deleted reference bases. Insertions
    mark their pileup anchor base so insertion evidence is read from ``PileupRead.indel``.
    """
    offsets = set()
    matcher = SequenceMatcher(a=wt_seq, b=design_donor, autojunk=False)
    for tag, i1, i2, _j1, _j2 in matcher.get_opcodes():
        if tag == "equal":
            continue
        if tag in ("replace", "delete"):
            offsets.update(range(i1, i2))
        elif tag == "insert" and wt_seq:
            offsets.add(min(max(i1 - 1, 0), len(wt_seq) - 1))
    return sorted(offsets)


# Max fraction of wt-vs-reference mismatches tolerated when deciding donor
# orientation. Real libraries carry strain SNPs + synthetic homology-arm recoding
# (observed ~5-13/129 = up to ~10% on the rearrayed NNS table, always reverse-
# complement), so an exact match is far too strict; random sequence sits near 0.75,
# so 0.35 cleanly separates a true (mis-SNP'd) match from genuinely wrong coords.
MAX_WT_REF_MISMATCH_FRAC = 0.35


def _hamming(a: str, b: str) -> int:
    return sum(1 for x, y in zip(a, b) if x != y)


def _site_design_offsets(fasta: "pysam.FastaFile", site: EditSite) -> Optional[List[int]]:
    """Designed offsets in genome orientation, or ``None`` for legacy whole-window scoring.

    Orientation (forward vs reverse-complement) is decided by the *better* match of
    ``wt_seq`` to the reference window, tolerating strain SNPs / synthetic recoding —
    not by exact equality. The designed offsets themselves come from the wt-vs-design
    diff and are independent of those reference mismatches. Raises ``ValueError`` when
    neither orientation matches the reference (wrong coords / wrong reference build).
    """
    if site.wt_seq is None or site.design_donor is None:
        return None

    ref = fasta.fetch(site.chrom, site.start - 1, site.end).upper()
    wt = site.wt_seq.upper()
    design = site.design_donor.upper()

    if len(wt) != len(ref):
        raise ValueError(
            f"{site.oligo or site.chrom}:{site.start}-{site.end} wt_seq length "
            f"{len(wt)} != reference window length {len(ref)}"
        )

    rc_wt = _rev_comp(wt)
    h_fwd = _hamming(wt, ref)
    h_rc = _hamming(rc_wt, ref)
    if min(h_fwd, h_rc) / len(ref) > MAX_WT_REF_MISMATCH_FRAC:
        raise ValueError(
            f"{site.oligo or site.chrom}:{site.start}-{site.end} wt_seq does not match "
            f"reference in either orientation (best {min(h_fwd, h_rc)}/{len(ref)} mismatches)"
        )

    if h_fwd <= h_rc:
        return _designed_offsets(wt, design)
    return _designed_offsets(rc_wt, _rev_comp(design))


def window_edit_vaf(
    bam: "pysam.AlignmentFile",
    fasta: "pysam.FastaFile",
    chrom: str,
    start: int,
    end: int,
    *,
    min_base_qual: int = 13,
) -> Tuple[float, int]:
    """Max non-reference allele fraction over the donor window ``[start, end]``
    (1-based inclusive), with the read depth at that position.

    This is the bridge-free edit signal: a real edit makes one or more positions in
    the window diverge from the reference at high fraction; the window maximum is
    robust to exactly where in the window the designed change sits. Deletions and
    ref-skips count as non-reference; insertions are scored at their anchor base.
    """
    ref = fasta.fetch(chrom, start - 1, end).upper()
    best_vaf, best_cov = 0.0, 0
    for col in bam.pileup(
        chrom, start - 1, end, truncate=True, stepper="samtools", min_base_quality=min_base_qual
    ):
        i = col.reference_pos - (start - 1)
        if i < 0 or i >= len(ref):
            continue
        refbase = ref[i]
        nonref = total = 0
        for pr in col.pileups:
            total += 1
            if pr.is_del or pr.is_refskip:
                nonref += 1
                continue
            qpos = pr.query_position
            if qpos is None:
                continue
            if pr.alignment.query_sequence[qpos].upper() != refbase or pr.indel > 0:
                nonref += 1
        if total == 0:
            continue
        vaf = nonref / total
        if vaf > best_vaf or (vaf == best_vaf and total > best_cov):
            best_vaf, best_cov = vaf, total
    return best_vaf, best_cov


def designed_site_vaf(
    bam: "pysam.AlignmentFile",
    fasta: "pysam.FastaFile",
    site: EditSite,
    *,
    min_base_qual: int = 13,
) -> Tuple[float, int]:
    """Edit VAF at the designed bases of ``site``.

    When ``wt_seq``/``design_donor`` are available, only designed offsets are
    scored. This is the low-background path used by the library-wide scan. Sites
    without sequence metadata fall back to the legacy donor-window maximum.
    """
    offsets = _site_design_offsets(fasta, site)
    if offsets is None:
        return window_edit_vaf(
            bam, fasta, site.chrom, site.start, site.end, min_base_qual=min_base_qual
        )
    if not offsets:
        return 0.0, 0

    wanted = set(offsets)
    ref = fasta.fetch(site.chrom, site.start - 1, site.end).upper()
    best_vaf, best_cov = 0.0, 0
    for col in bam.pileup(
        site.chrom,
        site.start - 1,
        site.end,
        truncate=True,
        stepper="samtools",
        min_base_quality=min_base_qual,
    ):
        i = col.reference_pos - (site.start - 1)
        if i not in wanted or i < 0 or i >= len(ref):
            continue
        refbase = ref[i]
        nonref = total = 0
        for pr in col.pileups:
            total += 1
            if pr.is_del or pr.is_refskip:
                nonref += 1
                continue
            qpos = pr.query_position
            if qpos is None:
                continue
            if pr.alignment.query_sequence[qpos].upper() != refbase or pr.indel > 0:
                nonref += 1
        if total == 0:
            continue
        vaf = nonref / total
        if vaf > best_vaf or (vaf == best_vaf and total > best_cov):
            best_vaf, best_cov = vaf, total
    return best_vaf, best_cov


def build_edit_sites(oligos: Iterable[str], oligo_windows: Dict[str, EditSite]) -> List[EditSite]:
    """Map bridged oligos to designed sites, de-duplicated by site identity."""
    seen: Dict[Tuple[str, int, int, Optional[str], Optional[str]], EditSite] = {}
    for oligo in oligos:
        site = oligo_windows.get(oligo)
        if site is None:
            continue
        seen.setdefault(site.key(), site)
    return list(seen.values())


def build_all_edit_sites(oligo_windows: Dict[str, EditSite]) -> List[EditSite]:
    """Return the library-wide designed-site set, de-duplicated across oligos."""
    seen: Dict[Tuple[str, int, int, Optional[str], Optional[str]], EditSite] = {}
    for site in oligo_windows.values():
        seen.setdefault(site.key(), site)
    return list(seen.values())


def load_oligo_windows(
    path: str,
    *,
    oligo_col: str = "oligo_name",
    chrom_col: str = "chrom",
    start_col: str = "donor_start_coord",
    end_col: str = "donor_end_coord",
    wt_col: str = "wt_seq",
    design_col: str = "synthetic_donor",
    design_col_fallback: str = "design_donor",
) -> Dict[str, EditSite]:
    """Load a per-library ``oligo_name -> EditSite`` map from a header'd TSV.

    The design column must be the **actually-synthesized / integrated** donor, since
    that is what appears in the edited genome. In the MAGESTIC tables this is
    ``synthetic_donor`` (target codon + synonymous re-cut-blocking mutations);
    ``design_donor`` records only the intended amino-acid change and is *identical to
    wt_seq for ~90% of oligos* (silent at the DNA level), so using it would make
    ``designed_site_vaf`` compute empty offsets and miss real edits. We therefore
    default to ``synthetic_donor`` and fall back to ``design_donor`` only if it is
    absent. When neither the wt nor design column is present, scoring falls back to
    the noisy donor-window maximum.
    """
    out: Dict[str, EditSite] = {}
    with open(path, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        idx = {name: i for i, name in enumerate(header)}
        oi, ci, si, ei = idx[oligo_col], idx[chrom_col], idx[start_col], idx[end_col]
        wi = idx.get(wt_col)
        di = idx.get(design_col)
        if di is None:
            di = idx.get(design_col_fallback)
        for row in reader:
            if len(row) <= max(oi, ci, si, ei):
                continue
            oligo = row[oi].strip()
            if not oligo or not row[si].strip().isdigit() or not row[ei].strip().isdigit():
                continue
            wt_seq = row[wi].strip().upper() if wi is not None and len(row) > wi and row[wi] else None
            design = row[di].strip().upper() if di is not None and len(row) > di and row[di] else None
            out.setdefault(
                oligo,
                EditSite(
                    row[ci].strip(),
                    int(row[si]),
                    int(row[ei]),
                    oligo,
                    wt_seq=wt_seq,
                    design_donor=design,
                ),
            )
    return out


# --------------------------------------------------------------------------- #
# Emit                                                                          #
# --------------------------------------------------------------------------- #

EMIT_COLUMNS = [
    "sample", "label", "confidence", "n_genome_bc", "n_plasmid_bc",
    "n_sites", "n_covered", "n_full", "n_partial", "n_wt", "n_nocov", "edits", "flags",
]


def reconciliation_to_row(c: CellReconciliation, *, emit_all_sites: bool = False) -> List[str]:
    detail_sites = c.sites if emit_all_sites else [
        s for s in c.sites if s.status in (SITE_FULL, SITE_PARTIAL)
    ]
    edits = ";".join(
        f"{s.site.oligo or '.'}@{s.site.chrom}:{s.site.start}-{s.site.end}={s.vaf:.2f}/{s.cov}:{s.status}"
        for s in detail_sites
    )
    return [
        c.sample, c.label, c.confidence, str(c.n_genome_bc), str(c.n_plasmid_bc),
        str(c.n_sites), str(c.n_covered), str(c.n_full), str(c.n_partial),
        str(c.n_wt), str(c.n_nocov), edits, ",".join(c.flags),
    ]


# --------------------------------------------------------------------------- #
# CLI                                                                           #
# --------------------------------------------------------------------------- #


def _build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="WGS Stage 3: reconcile bc1 lineage set x edit VAF -> "
        "SINGLE/MULTI/CONTAMINATION per cell."
    )
    p.add_argument("--bam-dir", required=True, help="directory of per-cell BAMs")
    p.add_argument("--bam-suffix", default=".sort.dedup.recal2.bam")
    p.add_argument("--fasta", required=True, help="reference fasta (indexed)")
    p.add_argument("--config", required=True, help="bc1 locus config (input.yaml subset)")
    p.add_argument("--bridge-table", required=True, help="bc1 -> oligo TSV (per library)")
    p.add_argument("--bridge-bc1-col", default="bc1")
    p.add_argument("--bridge-oligo-col", default="oligo_name")
    p.add_argument("--bridge-count-col", default=None)
    p.add_argument("--oligo-window-table", required=True, help="oligo -> donor window TSV")
    p.add_argument(
        "--design-col", default="synthetic_donor",
        help="oligo-window column holding the integrated donor sequence "
        "(synthetic_donor in MAGESTIC tables; design_donor is silent for ~90%% of oligos)",
    )
    p.add_argument("--wt-col", default="wt_seq", help="oligo-window column holding the wild-type seq")
    p.add_argument(
        "--site-scope",
        choices=("designed", "bridged"),
        default="designed",
        help=(
            "designed = scan every library designed site for bridge-independent edit counts; "
            "bridged = score only sites implied by confident bridged bc1s"
        ),
    )
    p.add_argument("--tau-full", type=float, default=0.85)
    p.add_argument("--tau-wt", type=float, default=0.15)
    p.add_argument("--min-cov", type=int, default=4)
    p.add_argument(
        "--min-partial-reads", type=int, default=3,
        help="absolute non-ref read-count floor for a PARTIAL (mixture) call; "
        "guards against single-read false PARTIALs across the designed-site scan",
    )
    p.add_argument("--min-support", type=int, default=2, help="min deduped molecules per bc1")
    p.add_argument(
        "--amplicon-left-window", default=None,
        help="genome-track amplicon-telltale left window, e.g. 10224-10235 (omit to disable)",
    )
    p.add_argument("--amplicon-right-window", default=None, help="e.g. 10318-10326")
    p.add_argument("--out", required=True)
    return p


def run(args: argparse.Namespace) -> int:
    if pysam is None:  # pragma: no cover
        raise RuntimeError("pysam is required for the BAM I/O layer")
    import glob
    import os

    from .bc1_consensus import (
        BridgeIndex,
        _parse_window,
        confident_barcodes,
        extract_track,
        load_locus_config,
        tracks_from_config,
    )

    cfg = load_locus_config(args.config)
    # The genome track needs its amplicon-telltale filter, or n_genome_bc is inflated
    # by the same amplicon/hop contamination that bc1_consensus filters out.
    genome_track, plasmid_track = tracks_from_config(
        cfg,
        amplicon_left_window=_parse_window(args.amplicon_left_window),
        amplicon_right_window=_parse_window(args.amplicon_right_window),
    )
    bridge = BridgeIndex.from_tsv(
        args.bridge_table, bc1_col=args.bridge_bc1_col,
        oligo_col=args.bridge_oligo_col, count_col=args.bridge_count_col,
    )
    oligo_windows = load_oligo_windows(
        args.oligo_window_table, wt_col=args.wt_col, design_col=args.design_col
    )

    bams = sorted(glob.glob(os.path.join(args.bam_dir, "*" + args.bam_suffix)))
    all_sites = build_all_edit_sites(oligo_windows)
    sys.stderr.write(
        f"[reconcile] {len(bams)} bams; {len(oligo_windows)} oligos; "
        f"{len(all_sites)} designed sites; site_scope={args.site_scope}\n"
    )
    if args.site_scope == "designed":
        # designed_site_vaf SILENTLY falls back to the noisy whole-window VAF
        # (~50% false-partial/site) when wt_seq/design_donor are missing; in the
        # library-wide designed scan that fallback would inflate false CONTAMINATION,
        # so surface it rather than letting it pass unnoticed.
        n_fallback = sum(1 for s in all_sites if s.wt_seq is None or s.design_donor is None)
        if n_fallback:
            sys.stderr.write(
                f"[reconcile] WARNING: {n_fallback}/{len(all_sites)} designed sites lack "
                "wt_seq/design_donor and will use the noisy whole-window VAF fallback\n"
            )

    with open(args.out, "w") as out, pysam.FastaFile(args.fasta) as fasta:
        out.write("\t".join(EMIT_COLUMNS) + "\n")

        # Pre-validate the designed-site set once against the reference (orientation is
        # BAM-independent): drop sites whose wt_seq matches neither orientation (wrong
        # coords / wrong reference build) with a visible count, rather than crashing the
        # whole run on the first one.
        designed_scan_sites: Optional[List[EditSite]] = None
        if args.site_scope == "designed":
            designed_scan_sites = []
            n_unmappable = 0
            for es in all_sites:
                try:
                    _site_design_offsets(fasta, es)
                except ValueError:
                    n_unmappable += 1
                    continue
                designed_scan_sites.append(es)
            if n_unmappable:
                sys.stderr.write(
                    f"[reconcile] WARNING: dropped {n_unmappable}/{len(all_sites)} designed "
                    "sites unmappable to the reference (wt_seq matches neither orientation)\n"
                )
            sys.stderr.write(
                f"[reconcile] scanning {len(designed_scan_sites)} reference-validated "
                "designed sites per cell\n"
            )

        for bampath in bams:
            sample = os.path.basename(bampath)[: -len(args.bam_suffix)]
            with pysam.AlignmentFile(bampath, "rb") as bam:
                gbcs = confident_barcodes(
                    extract_track(bam, genome_track), bridge,
                    min_support=args.min_support, high_only=True,
                )
                pbcs = (
                    confident_barcodes(
                        extract_track(bam, plasmid_track), bridge,
                        min_support=args.min_support, high_only=True,
                    )
                    if plasmid_track is not None
                    else []
                )
                # Edit count is bridge-independent in the default designed-site scan.
                # The old bridged scope is retained for narrow/debug runs.
                oligos = [b.oligo for b in (gbcs + pbcs) if b.oligo]
                edit_sites = (
                    designed_scan_sites
                    if designed_scan_sites is not None
                    else build_edit_sites(oligos, oligo_windows)
                )
                sites: List[SiteCall] = []
                for es in edit_sites:
                    try:
                        vaf, cov = designed_site_vaf(bam, fasta, es)
                    except ValueError:
                        # belt-and-suspenders for the bridged scope (designed scope is
                        # pre-validated above); skip a site that cannot be oriented.
                        continue
                    sites.append(
                        SiteCall(
                            es, vaf, cov,
                            classify_site(
                                vaf, cov, tau_full=args.tau_full,
                                tau_wt=args.tau_wt, min_cov=args.min_cov,
                                min_partial_reads=args.min_partial_reads,
                            ),
                        )
                    )
            # Confidence reflects the actual genome-track bc1 evidence in BOTH scopes.
            # In designed scope the *edit label* is bridge-independent (read from VAF at
            # designed bases), but the confidence field still tracks lineage/attribution
            # evidence, so it must not be forced to 'high' regardless of bc1 quality.
            bc_conf = "high" if (gbcs and gbcs[0].confidence in ("high", "medium")) else "low"
            cell = reconcile_cell(
                sample, sites, n_genome_bc=len(gbcs), n_plasmid_bc=len(pbcs), bc_confidence=bc_conf,
            )
            out.write("\t".join(reconciliation_to_row(cell)) + "\n")
    return 0


def main(argv: Optional[List[str]] = None) -> int:
    return run(_build_arg_parser().parse_args(argv))


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
