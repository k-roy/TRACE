"""Allele-specific WGS reconciliation and explicit window compatibility evidence.

All coordinates are 1-based inclusive in the BAM reference frame. Designed edits
require declared REF/ALT alleles; a non-reference maximum in a donor window is
compatibility evidence only. Intermediate intended-allele fractions and multiple
genome barcodes are mixture warnings, not a proof of biological contamination.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from numbers import Integral
import os
import tempfile
from contextlib import contextmanager
import sys
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

try:
    import pysam
except ImportError:  # pragma: no cover
    pysam = None

SITE_FULL = "full"
SITE_PARTIAL = "partial"
SITE_WT = "wt"
SITE_OTHER = "other"
SITE_NOCOV = "nocov"
SITE_UNASSESSED = "unassessed"
SITE_INCOMPLETE = "incomplete"
NO_EDIT = "NO_EDIT"
SINGLE_EDIT = "SINGLE_EDIT"
MULTI_EDIT = "MULTI_EDIT"
CONTAMINATION = "CONTAMINATION"
LOW_CONF = "LOW_CONF"
PARTIAL_EDIT = "PARTIAL_EDIT"


@dataclass(frozen=True)
class IntendedAllele:
    pos: int
    ref: str
    alt: str

    def __post_init__(self):
        if isinstance(self.pos, bool) or not isinstance(self.pos, int) or self.pos < 1:
            raise ValueError("intended allele pos must be a positive 1-based integer")
        if (not self.ref or not self.alt or set(self.ref) - set("ACGT")
                or set(self.alt) - set("ACGT") or self.ref == self.alt):
            raise ValueError("intended alleles require distinct, nonempty uppercase A/C/G/T REF and ALT")
        if len(self.ref) != len(self.alt) and self.ref[0] != self.alt[0]:
            raise ValueError("indels must include a shared left anchor base")


@dataclass(frozen=True)
class EditSite:
    chrom: str
    start: int
    end: int
    oligo: Optional[str] = None
    intended_alleles: Tuple[IntendedAllele, ...] = ()

    def __post_init__(self):
        if (not self.chrom or not isinstance(self.start, Integral) or not isinstance(self.end, Integral)
                or isinstance(self.start, bool) or isinstance(self.end, bool) or self.start < 1 or self.end < self.start):
            raise ValueError("edit window must have a contig and positive ordered coordinates")
        for allele in self.intended_alleles:
            if allele.pos < self.start or allele.pos + len(allele.ref) - 1 > self.end:
                raise ValueError("intended REF allele falls outside its donor window")
        spans = sorted((a.pos, a.pos + len(a.ref) - 1) for a in self.intended_alleles)
        if any(b[0] <= a[1] for a, b in zip(spans, spans[1:])):
            raise ValueError("overlapping intended alleles require one explicit combined allele")


@dataclass(frozen=True)
class AlleleEvidence:
    allele: IntendedAllele
    alt_count: int
    ref_count: int
    other_count: int
    depth: int
    low_mapq_reads: int = 0
    unknown_mapq_reads: int = 0

    @property
    def vaf(self):
        return self.alt_count / self.depth if self.depth else 0.0


@dataclass
class SiteCall:
    site: EditSite
    vaf: float
    cov: int
    status: str
    allele_evidence: Tuple[AlleleEvidence, ...] = ()
    flags: Tuple[str, ...] = ()


@dataclass
class CellReconciliation:
    sample: str
    label: str
    n_genome_bc: int
    n_plasmid_bc: int
    n_full: int
    n_partial: int
    n_nocov: int
    sites: List[SiteCall] = field(default_factory=list)
    confidence: str = "low"
    flags: Tuple[str, ...] = ()
    n_unassessed: int = 0
    n_incomplete: int = 0


def classify_site(vaf, cov, *, tau_full=0.85, tau_wt=0.15, min_cov=4):
    if (not (0 <= tau_wt < tau_full <= 1) or not isinstance(min_cov, Integral)
            or isinstance(min_cov, bool) or min_cov < 1):
        raise ValueError("require 0 <= tau_wt < tau_full <= 1 and min_cov >= 1")
    if (not math.isfinite(vaf) or not 0 <= vaf <= 1 or not isinstance(cov, Integral)
            or isinstance(cov, bool) or cov < 0):
        raise ValueError("VAF must be finite in [0, 1] and coverage nonnegative")
    if cov < min_cov:
        return SITE_NOCOV
    if vaf >= tau_full:
        return SITE_FULL
    if vaf <= tau_wt:
        return SITE_WT
    return SITE_PARTIAL


def bc_confidence_from(calls):
    """Barcode observation confidence is independent of bridge availability."""
    return "high" if calls and all(c.caller_confidence == "high" for c in calls) else "low"


def reconcile_cell(sample, sites, *, n_genome_bc=0, n_plasmid_bc=0, bc_confidence="low", flags=()):
    valid = {SITE_FULL, SITE_PARTIAL, SITE_WT, SITE_OTHER, SITE_NOCOV, SITE_UNASSESSED, SITE_INCOMPLETE}
    if any(s.status not in valid for s in sites):
        raise ValueError("unknown site status")
    n_full = sum(s.status == SITE_FULL for s in sites)
    n_partial = sum(s.status == SITE_PARTIAL for s in sites)
    n_nocov = sum(s.status == SITE_NOCOV for s in sites)
    n_unassessed = sum(s.status == SITE_UNASSESSED for s in sites)
    n_incomplete = sum(s.status == SITE_INCOMPLETE for s in sites)
    n_covered = sum(s.status in (SITE_FULL, SITE_PARTIAL, SITE_WT, SITE_OTHER, SITE_INCOMPLETE) for s in sites)
    result_flags = list(flags)
    result_flags.extend(flag for site in sites for flag in site.flags)
    if n_partial:
        result_flags.append("PARTIAL_VAF")
    if n_genome_bc >= 2:
        result_flags.append("MULTI_GENOME_BC")
    if any(s.status == SITE_OTHER or any(e.other_count for e in s.allele_evidence) for s in sites):
        result_flags.append("OTHER_ALLELE")
    if n_nocov:
        result_flags.append("UNCOVERED_SITES")
    if n_unassessed:
        result_flags.append("INTENDED_ALLELES_UNSPECIFIED")
    if n_incomplete:
        result_flags.append("INCOMPLETE_DESIGN")
    if n_partial or n_genome_bc >= 2:
        label = CONTAMINATION
        result_flags.append("MIXTURE_SIGNAL_NOT_PROOF")
    elif n_incomplete:
        label = PARTIAL_EDIT
    elif n_full >= 2:
        label = MULTI_EDIT
    elif n_full:
        label = SINGLE_EDIT
    elif n_covered:
        label = NO_EDIT
    else:
        label = LOW_CONF
    confidence = "high" if bc_confidence == "high" else "low"
    if label == LOW_CONF or n_nocov or n_unassessed or n_incomplete or "UNRESOLVED_BC1" in result_flags or "UNCOVERED_ALLELES" in result_flags:
        confidence = "low"
    return CellReconciliation(sample, label, n_genome_bc, n_plasmid_bc, n_full, n_partial,
                              n_nocov, sites, confidence, tuple(dict.fromkeys(result_flags)),
                              n_unassessed, n_incomplete)


def window_edit_vaf(bam, fasta, chrom, start, end, *, min_base_qual=13, min_cov=4, min_mapq=20):
    """Window compatibility maximum, preferring positions with sufficient coverage.

    This legacy observable does not establish an intended edit. Ref-skips are
    missing observations, never alternate alleles. If no position reaches min_cov,
    return the deepest observable position so low coverage remains explicit.
    """
    ref = fasta.fetch(chrom, start - 1, end).upper()
    best = (0.0, 0)
    covered = False
    for col in bam.pileup(chrom, start - 1, end, truncate=True, stepper="samtools",
                          min_base_quality=min_base_qual, flag_filter=0xF04):
        i = col.reference_pos - start + 1
        if i < 0 or i >= len(ref):
            continue
        nonref = total = 0
        for pr in col.pileups:
            if pr.alignment.mapping_quality == 255 or pr.alignment.mapping_quality < min_mapq:
                continue
            if pr.is_refskip:
                continue
            if pr.is_del:
                total += 1
                nonref += 1
                continue
            qpos = pr.query_position
            if qpos is None:
                continue
            base = pr.alignment.query_sequence[qpos].upper()
            if base not in "ACGT":
                continue
            total += 1
            nonref += base != ref[i] or pr.indel > 0
        if not total:
            continue
        vaf = nonref / total
        if total >= min_cov:
            if not covered or (vaf, total) > best:
                best = (vaf, total)
            covered = True
        elif not covered and (total, vaf) > (best[1], best[0]):
            best = (vaf, total)
    return best


def _read_allele(read, allele, min_base_qual):
    """Exact read haplotype across REF, including insertions after REF bases.

    Return None for incomplete span, skipped reference, ambiguous/low-quality
    bases, or a deletion of the anchor. Hard clips provide no sequence.
    """
    start, end = allele.pos - 1, allele.pos - 1 + len(allele.ref)
    if read.reference_start > start or read.reference_end is None or read.reference_end < end:
        return None
    seq, qualities = read.query_sequence, read.query_qualities
    if seq is None or qualities is None:
        return None
    q, r = 0, read.reference_start
    observed = []
    anchor_seen = False
    need_right_anchor = len(allele.ref) != len(allele.alt)
    right_seen = not need_right_anchor
    for op, size in read.cigartuples or ():
        if op in (0, 7, 8):
            if need_right_anchor and r <= end < r + size:
                right_q = q + end - r
                right_seen = seq[right_q].upper() in "ACGT" and qualities[right_q] >= min_base_qual
            lo, hi = max(r, start), min(r + size, end)
            if lo < hi:
                qlo, qhi = q + lo - r, q + hi - r
                bases = seq[qlo:qhi].upper()
                if set(bases) - set("ACGT") or any(x < min_base_qual for x in qualities[qlo:qhi]):
                    return None
                observed.append(bases)
                anchor_seen |= lo == start
            q += size
            r += size
        elif op == 1:
            if start < r <= end:
                bases = seq[q:q + size].upper()
                if set(bases) - set("ACGT") or any(x < min_base_qual for x in qualities[q:q + size]):
                    return None
                observed.append(bases)
            q += size
        elif op in (2, 3):
            if r < end and r + size > start:
                if op == 3 or r <= start or r + size > end:
                    return None
            r += size
        elif op == 4:
            q += size
        elif op not in (5, 6):
            return None
        if r > end:
            break
    return "".join(observed) if anchor_seen and right_seen else None


def _normalize_allele(fasta, chrom, allele):
    """Minimal left-aligned representation, using bounded reference chunks."""
    pos, ref, alt = allele.pos, allele.ref, allele.alt
    cached_start, cached = -1, ""
    while ref[-1] == alt[-1]:
        if min(len(ref), len(alt)) == 1:
            if pos == 1:
                break
            offset = pos - 2
            if not cached_start <= offset < cached_start + len(cached):
                cached_start = max(0, offset - 4095)
                cached = fasta.fetch(chrom, cached_start, offset + 1).upper()
            previous = cached[offset - cached_start]
            if previous not in "ACGT":
                break
            ref, alt, pos = previous + ref, previous + alt, pos - 1
        ref, alt = ref[:-1], alt[:-1]
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    return (pos, ref, alt)


def _equivalent_indel(read, fasta, chrom, allele, normalized, min_base_qual, normalizations):
    """Validate a reference-equivalent indel against the observed local haplotype.

    Reference normalization identifies candidate events; it cannot replace an
    observed mutant anchor or move an indel across an observed repeat-breaking
    SNP. Compare the complete interval joining the requested and candidate events
    to the reference haplotype with exactly the intended allele applied. Include
    every equivalent candidate in that interval so two indels cannot masquerade
    as one. Return alt/other/unknown, or None when no event needs this check.
    """
    seq, qualities = read.query_sequence, read.query_qualities
    if seq is None or qualities is None:
        return None
    q, pos = 0, read.reference_start
    # Padding consumes neither coordinate; genomic neighbors remain adjacent.
    cigar = [(op, size) for op, size in read.cigartuples or () if op != 6]
    candidates = []
    for index, (op, size) in enumerate(cigar):
        if op in (0, 7, 8):
            q += size
            pos += size
        elif op in (1, 2):
            right = q + size if op == 1 else q
            aligned_neighbors = (index > 0 and index + 1 < len(cigar)
                                 and cigar[index - 1][0] in (0, 7, 8)
                                 and cigar[index + 1][0] in (0, 7, 8))
            if aligned_neighbors and pos > 0 and q > 0 and right < len(seq):
                checked = seq[q - 1:right + 1].upper()
                if not set(checked) - set("ACGT") and min(qualities[q - 1:right + 1]) >= min_base_qual:
                    ref = fasta.fetch(chrom, pos - 1, pos + (size if op == 2 else 0)).upper()
                    alt = ref[0] + (seq[q:q + size].upper() if op == 1 else "")
                    if not set(ref) - set("ACGT"):
                        key = (pos, ref, alt)
                        if key not in normalizations:
                            normalizations[key] = _normalize_allele(fasta, chrom, IntendedAllele(*key))
                        if normalizations[key] == normalized:
                            candidates.append((pos, pos + len(ref) - 1))
            if op == 1:
                q += size
            else:
                pos += size
        elif op == 3:
            pos += size
        elif op == 4:
            q += size
    if not candidates:
        return None
    start = min(allele.pos, min(a for a, b in candidates))
    end = max(allele.pos + len(allele.ref) - 1, max(b for a, b in candidates))
    reference = fasta.fetch(chrom, start - 1, end).upper()
    offset = allele.pos - start
    intended = reference[:offset] + allele.alt + reference[offset + len(allele.ref):]
    observation = _read_allele(read, IntendedAllele(start, reference, intended), min_base_qual)
    if observation is None:
        return "unknown"
    return "alt" if observation == intended else "other"


def _trim_allele_padding(allele):
    """Remove unchanged representation padding without shifting indels in repeats.

    At least one REF/ALT base remains, retaining the required indel anchor. This
    limits the measurement interval to the declared change; left normalization
    remains a separate equivalence operation and never moves the BAM fetch away
    from the observed locus merely because its reference is repetitive.
    """
    pos, ref, alt = allele.pos, allele.ref, allele.alt
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        # A complex delins also needs its final shared left anchor. Pure
        # substitutions do not, so remove all their unchanged prefix bases.
        if len(ref) != len(alt) and ref[1] != alt[1]:
            break
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    return IntendedAllele(pos, ref, alt)


def intended_allele_evidence(bam, fasta, chrom, allele, *, min_base_qual=13, min_mapq=20):
    """Count independent read fragments supporting exact intended REF/ALT.

    Primary, mapped, nonduplicate, QC-passing reads with MAPQ >= min_mapq are used.
    MAPQ 255 (unknown) is always excluded. Rejected overlapping primary alignment
    records are counted separately as low_mapq_reads and unknown_mapq_reads. Overlapping mates
    with one read name/read group count once; discordant mate observations are
    discarded. Unchanged REF/ALT padding is removed for measurement; the original
    declaration is retained in the returned evidence. Counts are fragments, unlike
    the legacy window pileup depth.
    """
    observed_ref = fasta.fetch(chrom, allele.pos - 1, allele.pos - 1 + len(allele.ref)).upper()
    if observed_ref != allele.ref:
        raise ValueError(f"{chrom}:{allele.pos}: intended REF {allele.ref!r} does not match reference {observed_ref!r}")
    declared_allele = allele
    allele = _trim_allele_padding(allele)
    fragments = {}
    normalized = _normalize_allele(fasta, chrom, allele) if len(allele.ref) != len(allele.alt) else None
    normalizations = {}
    low_mapq_reads = unknown_mapq_reads = 0
    for index, read in enumerate(bam.fetch(chrom, allele.pos - 1, allele.pos - 1 + len(allele.ref))):
        if read.flag & 0xF04:
            continue
        if read.mapping_quality == 255:
            unknown_mapq_reads += 1
            continue
        if read.mapping_quality < min_mapq:
            low_mapq_reads += 1
            continue
        obs = _read_allele(read, allele, min_base_qual)
        if (normalized is not None and min(len(normalized[1]), len(normalized[2])) == 1
                and read.reference_start <= allele.pos - 1
                and read.reference_end is not None
                and read.reference_end >= allele.pos - 1 + len(allele.ref)):
            equivalent = _equivalent_indel(read, fasta, chrom, allele, normalized,
                                           min_base_qual, normalizations)
            if equivalent == "alt":
                obs = allele.alt
            elif equivalent == "other":
                obs = "OTHER_INDEL_HAPLOTYPE"
            elif equivalent == "unknown":
                obs = None
        if obs is None:
            continue
        # Unpaired reads are distinct observations even when their names coincide.
        key = (read.get_tag("RG") if read.has_tag("RG") else "", read.query_name) if read.is_paired else (None, index)
        if key in fragments and fragments[key] != obs:
            fragments[key] = None
        else:
            fragments[key] = obs
    values = [x for x in fragments.values() if x is not None]
    alt = sum(x == allele.alt for x in values)
    ref = sum(x == allele.ref for x in values)
    return AlleleEvidence(declared_allele, alt, ref, len(values) - alt - ref, len(values), low_mapq_reads, unknown_mapq_reads)


def score_edit_site(bam, fasta, site, *, tau_full=0.85, tau_wt=0.15, min_cov=4, min_base_qual=13, min_mapq=20):
    classify_site(0, 0, tau_full=tau_full, tau_wt=tau_wt, min_cov=min_cov)
    if not isinstance(min_mapq, Integral) or isinstance(min_mapq, bool) or not 0 <= min_mapq <= 254:
        raise ValueError("minimum mapping quality must be an integer from 0 to 254")
    if not isinstance(min_base_qual, Integral) or min_base_qual < 0:
        raise ValueError("minimum base quality must be a nonnegative integer")
    if not site.intended_alleles:
        vaf, cov = window_edit_vaf(bam, fasta, site.chrom, site.start, site.end,
                                   min_base_qual=min_base_qual, min_cov=min_cov, min_mapq=min_mapq)
        return SiteCall(site, vaf, cov, SITE_UNASSESSED)
    evidence = tuple(intended_allele_evidence(bam, fasta, site.chrom, allele,
                                             min_base_qual=min_base_qual, min_mapq=min_mapq) for allele in site.intended_alleles)
    statuses = [classify_site(e.vaf, e.depth, tau_full=tau_full, tau_wt=tau_wt, min_cov=min_cov) for e in evidence]
    if SITE_PARTIAL in statuses:
        status = SITE_PARTIAL
    elif SITE_NOCOV in statuses:
        status = SITE_NOCOV
    elif all(s == SITE_FULL for s in statuses):
        status = SITE_FULL
    elif all(s == SITE_WT for s in statuses):
        status = SITE_OTHER if any(e.other_count / e.depth > tau_wt for e in evidence) else SITE_WT
    else:
        status = SITE_INCOMPLETE
    # Minimum intended-allele fraction/depth is an explicit conservative summary,
    # not a phased all-edits-on-one-molecule fraction. Exact counts are also emitted.
    flags = []
    if SITE_NOCOV in statuses:
        flags.append("UNCOVERED_ALLELES")
    if any(e.low_mapq_reads or e.unknown_mapq_reads for e in evidence):
        flags.append("AMBIGUOUS_ALIGNMENTS_EXCLUDED")
    return SiteCall(site, min(e.vaf for e in evidence), min(e.depth for e in evidence), status, evidence, tuple(flags))


def build_edit_sites(oligos, oligo_windows, *, fasta=None):
    """One site per distinct intended event set; reject overlapping hypotheses.

    With a reference, repeat-equivalent and padded allele representations are
    normalized before deciding whether two donors describe the same event.
    """
    seen = {}
    for oligo in oligos:
        if oligo not in oligo_windows:
            raise ValueError(f"no donor window for assigned/bridged oligo {oligo!r}")
        win = oligo_windows[oligo]
        site = win if isinstance(win, EditSite) else EditSite(*win, oligo)
        identity = []
        for allele in site.intended_alleles:
            if fasta is not None:
                observed = fasta.fetch(site.chrom, allele.pos - 1, allele.pos - 1 + len(allele.ref)).upper()
                if observed != allele.ref:
                    raise ValueError(f"{site.chrom}:{allele.pos}: intended REF does not match reference")
                identity.append(_normalize_allele(fasta, site.chrom, allele))
            else:
                identity.append((allele.pos, allele.ref, allele.alt))
        identity = tuple(sorted(identity))
        key = ((site.chrom, identity) if identity else (site.chrom, site.start, site.end))
        for previous_key, previous in seen.items():
            if previous.chrom != site.chrom:
                continue
            previous_identity = previous_key[1] if len(previous_key) == 2 else ()
            if identity and previous_identity:
                if identity == previous_identity:
                    continue
                spans = [(pos, pos + len(ref) - 1) for pos, ref, alt in previous_identity]
                if any(pos <= right and pos + len(ref) - 1 >= left
                       for pos, ref, alt in identity for left, right in spans):
                    raise ValueError(f"overlapping intended designs at {site.chrom}; cannot count independent edits")
            if ((site.start, site.end) == (previous.start, previous.end)
                    and identity != previous_identity):
                raise ValueError(f"conflicting intended designs at {site.chrom}:{site.start}-{site.end}")
        seen.setdefault(key, site)
    return list(seen.values())


def load_oligo_windows(path, *, oligo_col="oligo_name", chrom_col="chrom",
                       start_col="donor_start_coord", end_col="donor_end_coord",
                       intended_col="intended_variants"):
    """Strict TSV loader; optional intended_variants is JSON [{pos,ref,alt}, ...].

    Alleles must be explicit intended changes in this reference coordinate frame.
    Donor/reference differences alone are not equivalent to intended edits.
    """
    out = {}
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        required = {oligo_col, chrom_col, start_col, end_col}
        if not reader.fieldnames or not required.issubset(reader.fieldnames) or len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise ValueError(f"{path}: missing or duplicate window header; require {sorted(required)}")
        for lineno, row in enumerate(reader, 2):
            if None in row or any(row.get(c) is None for c in required):
                raise ValueError(f"{path}:{lineno}: malformed donor-window row")
            oligo = row[oligo_col].strip()
            try:
                if not oligo:
                    raise ValueError("empty oligo name")
                raw = row.get(intended_col, "") or ""
                variants = json.loads(raw) if raw.strip() else []
                if not isinstance(variants, list) or any(not isinstance(v, dict) or set(v) != {"pos", "ref", "alt"} for v in variants):
                    raise ValueError("intended_variants must be a list of {pos,ref,alt} objects")
                alleles = tuple(sorted((IntendedAllele(**v) for v in variants), key=lambda a: (a.pos, a.ref, a.alt)))
                site = EditSite(row[chrom_col].strip(), int(row[start_col]), int(row[end_col]), oligo, alleles)
            except (TypeError, ValueError) as exc:
                raise ValueError(f"{path}:{lineno}: invalid window: {exc}") from exc
            if oligo in out and out[oligo] != site:
                raise ValueError(f"{path}:{lineno}: conflicting donor windows for {oligo!r}")
            out[oligo] = site
    if not out:
        raise ValueError(f"{path}: donor-window table has no data rows")
    return out


EMIT_COLUMNS = ["sample", "label", "confidence", "n_genome_bc", "n_plasmid_bc", "n_full",
                "n_partial", "n_nocov", "edits", "flags", "n_unassessed", "n_incomplete", "allele_counts"]


def reconciliation_to_row(c):
    edits = ";".join(f"{s.site.oligo or '.'}@{s.site.chrom}:{s.site.start}-{s.site.end}={s.vaf:.2f}/{s.cov}:{s.status}" for s in c.sites)
    counts = [{"oligo": s.site.oligo, "chrom": s.site.chrom, "pos": e.allele.pos,
               "ref": e.allele.ref, "alt": e.allele.alt, "alt_count": e.alt_count,
               "ref_count": e.ref_count, "other_count": e.other_count, "depth": e.depth,
               "low_mapq_reads": e.low_mapq_reads, "unknown_mapq_reads": e.unknown_mapq_reads}
              for s in c.sites for e in s.allele_evidence]
    return [c.sample, c.label, c.confidence, str(c.n_genome_bc), str(c.n_plasmid_bc), str(c.n_full),
            str(c.n_partial), str(c.n_nocov), edits, ",".join(c.flags), str(c.n_unassessed),
            str(c.n_incomplete), json.dumps(counts, separators=(",", ":"))]


@contextmanager
def _atomic_output(path):
    """An exception preserves any old complete output and removes the new prefix."""
    directory = os.path.dirname(os.path.abspath(path))
    fd, temporary = tempfile.mkstemp(prefix=".reconcile-", dir=directory, text=True)
    try:
        with os.fdopen(fd, "w") as out:
            yield out
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def _build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="WGS Stage 3: reconcile bc1 lineage set x edit VAF -> "
        "SINGLE/MULTI/CONTAMINATION per cell."
    )
    p.add_argument("--bam-dir", required=True, help="directory of per-cell BAMs")
    p.add_argument("--bam-suffix", default=".bam")
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
    p.add_argument("--min-mapq", type=int, default=20, help="variant evidence floor; 255 always excluded")
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

    classify_site(0, 0, tau_full=args.tau_full, tau_wt=args.tau_wt, min_cov=args.min_cov)
    if not args.bam_suffix:
        raise ValueError("BAM suffix must be nonempty")
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
    if not bams:
        raise ValueError(f"no BAMs match {args.bam_dir!r} and suffix {args.bam_suffix!r}")
    sys.stderr.write(f"[reconcile] {len(bams)} bams; {len(oligo_windows)} oligo windows\n")

    with _atomic_output(args.out) as out, pysam.FastaFile(args.fasta) as fasta:
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
                for es in build_edit_sites(oligos, oligo_windows, fasta=fasta):
                    sites.append(score_edit_site(bam, fasta, es, tau_full=args.tau_full,
                                                 tau_wt=args.tau_wt, min_cov=args.min_cov, min_mapq=args.min_mapq))
            bc_conf = bc_confidence_from(gbcs)
            cell = reconcile_cell(
                sample, sites, n_genome_bc=len(gbcs), n_plasmid_bc=len(pbcs), bc_confidence=bc_conf,
                flags=("UNRESOLVED_BC1",) if any(not b.oligo for b in gbcs + pbcs) else (),
            )
            out.write("\t".join(reconciliation_to_row(cell)) + "\n")
    return 0


def main(argv: Optional[List[str]] = None) -> int:
    return run(_build_arg_parser().parse_args(argv))


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
