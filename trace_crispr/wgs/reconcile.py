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
from genomic VAF at designed sites** (robust, bridge-free), while bc1 bridging is
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
from typing import Dict, List, Optional, Tuple

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
    """

    chrom: str
    start: int
    end: int
    oligo: Optional[str] = None


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
    n_full: int
    n_partial: int
    n_nocov: int
    sites: List[SiteCall] = field(default_factory=list)
    confidence: str = "low"  # high / low
    flags: Tuple[str, ...] = ()


# --------------------------------------------------------------------------- #
# Pure classification core (no BAM — fully unit-testable)                      #
# --------------------------------------------------------------------------- #


def classify_site(
    vaf: float, cov: int, *, tau_full: float = 0.85, tau_wt: float = 0.15, min_cov: int = 4
) -> str:
    """Bucket one designed site by its edit VAF and coverage."""
    if cov < min_cov:
        return SITE_NOCOV
    if vaf >= tau_full:
        return SITE_FULL
    if vaf <= tau_wt:
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
    n_full = sum(1 for s in sites if s.status == SITE_FULL)
    n_partial = sum(1 for s in sites if s.status == SITE_PARTIAL)
    n_nocov = sum(1 for s in sites if s.status == SITE_NOCOV)
    n_covered = sum(1 for s in sites if s.status != SITE_NOCOV)

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
        n_full=n_full, n_partial=n_partial, n_nocov=n_nocov, sites=sites,
        confidence=confidence, flags=tuple(flags),
    )


# --------------------------------------------------------------------------- #
# BAM I/O layer                                                                #
# --------------------------------------------------------------------------- #


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


def build_edit_sites(
    oligos: List[str], oligo_windows: Dict[str, Tuple[str, int, int]]
) -> List[EditSite]:
    """Map a cell's bridged oligos to their designed donor windows, de-duplicated by
    genomic window so two bc1 bridging to the same oligo/site are scored once."""
    seen: Dict[Tuple[str, int, int], EditSite] = {}
    for oligo in oligos:
        win = oligo_windows.get(oligo)
        if win is None:
            continue
        chrom, start, end = win
        seen.setdefault((chrom, start, end), EditSite(chrom, start, end, oligo))
    return list(seen.values())


def load_oligo_windows(
    path: str,
    *,
    oligo_col: str = "oligo_name",
    chrom_col: str = "chrom",
    start_col: str = "donor_start_coord",
    end_col: str = "donor_end_coord",
) -> Dict[str, Tuple[str, int, int]]:
    """Load a per-library ``oligo_name -> (chrom, start, end)`` map from a header'd
    TSV (e.g. Shengdi's ``sample_BC_variant_fil`` or an oligo design table)."""
    out: Dict[str, Tuple[str, int, int]] = {}
    with open(path, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        idx = {name: i for i, name in enumerate(header)}
        oi, ci, si, ei = idx[oligo_col], idx[chrom_col], idx[start_col], idx[end_col]
        for row in reader:
            if len(row) <= max(oi, ci, si, ei):
                continue
            oligo = row[oi].strip()
            if not oligo or not row[si].strip().isdigit() or not row[ei].strip().isdigit():
                continue
            out.setdefault(oligo, (row[ci].strip(), int(row[si]), int(row[ei])))
    return out


# --------------------------------------------------------------------------- #
# Emit                                                                          #
# --------------------------------------------------------------------------- #

EMIT_COLUMNS = [
    "sample", "label", "confidence", "n_genome_bc", "n_plasmid_bc",
    "n_full", "n_partial", "n_nocov", "edits", "flags",
]


def reconciliation_to_row(c: CellReconciliation) -> List[str]:
    edits = ";".join(
        f"{s.site.oligo or '.'}@{s.site.chrom}:{s.site.start}-{s.site.end}={s.vaf:.2f}/{s.cov}:{s.status}"
        for s in c.sites
    )
    return [
        c.sample, c.label, c.confidence, str(c.n_genome_bc), str(c.n_plasmid_bc),
        str(c.n_full), str(c.n_partial), str(c.n_nocov), edits, ",".join(c.flags),
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
    p.add_argument("--tau-full", type=float, default=0.85)
    p.add_argument("--tau-wt", type=float, default=0.15)
    p.add_argument("--min-cov", type=int, default=4)
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
    oligo_windows = load_oligo_windows(args.oligo_window_table)

    bams = sorted(glob.glob(os.path.join(args.bam_dir, "*" + args.bam_suffix)))
    sys.stderr.write(f"[reconcile] {len(bams)} bams; {len(oligo_windows)} oligo windows\n")

    with open(args.out, "w") as out, pysam.FastaFile(args.fasta) as fasta:
        out.write("\t".join(EMIT_COLUMNS) + "\n")
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
                # Edit set = the union of all confident (genome + plasmid) bridged
                # oligos' designed windows; the EDIT COUNT then comes from VAF, not
                # from how many bc1 bridged.
                oligos = [b.oligo for b in (gbcs + pbcs) if b.oligo]
                sites: List[SiteCall] = []
                for es in build_edit_sites(oligos, oligo_windows):
                    vaf, cov = window_edit_vaf(bam, fasta, es.chrom, es.start, es.end)
                    sites.append(
                        SiteCall(
                            es, vaf, cov,
                            classify_site(
                                vaf, cov, tau_full=args.tau_full,
                                tau_wt=args.tau_wt, min_cov=args.min_cov,
                            ),
                        )
                    )
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
