"""
trace_crispr.wgs.donor_outcome
==============================

On-target editing-outcome **pre-calling** for WGS MAGESTIC data.

For each GATK-called variant inside a sample's donor window, decide whether it is
a *designed* edit, a *synthetic error* (introduced during oligo synthesis /
cloning — i.e. present in the step-2 "synthetic" donor), or a *spontaneous*
mutation, using the four-Levenshtein bucketing of Shengdi Li, then aggregate the
per-variant calls into a per-sample editing outcome
(``designed_variant_installed`` / ``..._partially_installed`` / ``..._unedited``
/ ``unknown``, with optional ``;synthetic_error_installed`` /
``;spontaneous_mutations`` overlays).

Four distances for a candidate variant (``donor_tmp`` = ``donor_wt`` with the
variant substituted in)::

    e1 = lev(donor_tmp, donor_design)    e2 = lev(donor_wt, donor_design)
    e3 = lev(donor_tmp, donor_step2)     e4 = lev(donor_wt, donor_step2)

    e1 <  e2               -> designed_variant     (moves toward the design)
    e1 >= e2 and e3 <  e4  -> synthetic_error      (toward the step-2 donor)
    e1 >= e2 and e3 >= e4  -> spontaneous_variant

This is the step2/synthetic axis that TRACE's amplicon classifier does not have;
the same :func:`bucket_variant` utility is reusable from amplicon classification.

Provenance
----------
Original: Shengdi Li, ``precalling_target_outcome.pl`` (2023-10-31), in
``WGS_MAGESTIC_QTL_outcome_mapping`` — a light version of the precalling step of
the MAGESTIC-SCORE workflow (https://github.com/shli-embl/MAGESTIC-SCORE).
Python port: Kevin R. Roy (2026), validated for bit-for-bit parity against
Shengdi's reference ``on_target_precalling.txt`` outputs.

The bucketing comparisons and omission of spontaneous low-confidence calls retain
legacy semantics. Parsing and safety intentionally differ: FORMAT is interpreted by
name, genotype alleles are selected explicitly, reference blocks require MIN_DP,
coordinates are inclusive, variants remain structured, missing inputs fail, and final
output is atomic. FILTER remains unfiltered for compatibility; consumers must choose
an appropriate called-variant input upstream. These contract repairs do not establish
biological accuracy of the four-distance decision rule.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import os
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

try:
    from rapidfuzz.distance import Levenshtein

    def levenshtein(s1: str, s2: str) -> int:
        """Unit-cost Levenshtein edit distance (== Shengdi's hand-rolled Perl impl)."""
        return Levenshtein.distance(s1, s2)

except ImportError:  # pragma: no cover - rapidfuzz is a trace-crispr dependency

    def levenshtein(s1: str, s2: str) -> int:
        """Pure-Python fallback unit-cost Levenshtein (used only without rapidfuzz)."""
        if not s1:
            return len(s2)
        if not s2:
            return len(s1)
        prev = list(range(len(s2) + 1))
        for i, c1 in enumerate(s1, 1):
            cur = [i]
            for j, c2 in enumerate(s2, 1):
                cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (c1 != c2)))
            prev = cur
        return prev[-1]


_COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N", "-": "-"}


def reverse_complement(seq: str) -> str:
    """Reverse complement of an uppercase DNA string (mirrors the Perl helper)."""
    out = []
    for base in reversed(seq):
        if base not in _COMPLEMENT:
            raise ValueError(f"non-recognized nucleotide letter {base!r}")
        out.append(_COMPLEMENT[base])
    return "".join(out)


def load_fasta(path: str) -> Dict[str, str]:
    """Load a (small) reference FASTA into ``{display_id: uppercase_sequence}``.

    ``display_id`` is the first whitespace-delimited token of the header, matching
    BioPerl's ``$seq->display_id``.
    """
    refs: Dict[str, str] = {}
    name: Optional[str] = None
    chunks: List[str] = []
    def add_contig():
        sequence = "".join(chunks).upper()
        if name in refs or not sequence:
            raise ValueError(f"duplicate or empty FASTA contig {name!r}")
        if set(sequence) - set("ACGTRYSWKMBDHVN"):
            raise ValueError(f"invalid DNA sequence in FASTA contig {name!r}")
        refs[name] = sequence
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    add_contig()
                if not line[1:].split():
                    raise ValueError("empty FASTA identifier")
                name = line[1:].split()[0]
                chunks = []
            else:
                if name is None and line.strip():
                    raise ValueError(f"sequence before first FASTA header: {path}")
                chunks.append(line.strip())
    if name is not None:
        add_contig()
    if not refs:
        raise ValueError(f"empty reference FASTA: {path}")
    return refs


def _reference_identity(path: Path) -> tuple:
    stat = path.stat()
    return (stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns)


@dataclass(frozen=True)
class LoadedReference:
    """One immutable-source FASTA load for a cohort sharing that reference."""

    path: Path
    identity: tuple
    sequences: Dict[str, str]

    def require_current(self, path: Path) -> Dict[str, str]:
        if path.resolve() != self.path or _reference_identity(path) != self.identity:
            raise ValueError("preloaded reference does not match current genome source")
        return self.sequences


def load_reference(path: Path) -> LoadedReference:
    path = Path(path).resolve()
    identity = _reference_identity(path)
    sequences = load_fasta(str(path))
    if identity != _reference_identity(path):
        raise ValueError("reference changed while loading")
    return LoadedReference(path, identity, sequences)


# --------------------------------------------------------------------------- #
# Core: four-Levenshtein bucketing (the reusable utility)                      #
# --------------------------------------------------------------------------- #

# Variant bucket labels
DESIGNED = "designed"
SYNTHETIC_ERROR = "synthetic_error"
SPONTANEOUS = "spontaneous"


def bucket_variant(
    donor_tmp: str,
    donor_wt: str,
    donor_design: str,
    donor_step2: str,
) -> Tuple[str, int, int]:
    """Classify a single candidate donor allele via the four-Levenshtein test.

    Args:
        donor_tmp:    ``donor_wt`` with the candidate variant applied.
        donor_wt:     wild-type donor-window sequence (genome orientation).
        donor_design: intended (designed) donor sequence.
        donor_step2:  step-2 / synthetic donor sequence.

    Returns:
        ``(bucket, e2, e1)`` where ``bucket`` is one of ``DESIGNED`` /
        ``SYNTHETIC_ERROR`` / ``SPONTANEOUS`` and ``e2 -> e1`` is the
        wt-vs-tmp distance-to-design pair Shengdi records as the call evidence.
    """
    e1 = levenshtein(donor_tmp, donor_design)
    e2 = levenshtein(donor_wt, donor_design)
    if e1 < e2:
        return DESIGNED, e2, e1
    e3 = levenshtein(donor_tmp, donor_step2)
    e4 = levenshtein(donor_wt, donor_step2)
    if e3 < e4:
        return SYNTHETIC_ERROR, e2, e1
    return SPONTANEOUS, e2, e1


def apply_variant(donor_wt: str, pos0: int, ref: str, alt: str) -> str:
    """Return ``donor_wt`` with ``ref`` at 0-based ``pos0`` replaced by ``alt``.

    Mirrors Perl ``substr($donor_tmp, $pos0, length($ref)) = $alt``.
    """
    if pos0 < 0 or pos0 + len(ref) > len(donor_wt) or not ref:
        raise ValueError(f"variant outside donor window: {pos0}, {ref!r}")
    if donor_wt[pos0:pos0 + len(ref)] != ref:
        raise ValueError(f"variant REF {ref!r} does not match donor at offset {pos0}")
    if not alt or set(ref + alt) - set("ACGTN"):
        raise ValueError(f"unsupported non-sequence allele: {ref!r} -> {alt!r}")
    return donor_wt[:pos0] + alt + donor_wt[pos0 + len(ref) :]


# --------------------------------------------------------------------------- #
# Per-sample pre-calling                                                       #
# --------------------------------------------------------------------------- #


@dataclass
class DonorSpec:
    """A sample's donor-window specification (one samplesheet row).

    donor_strand declares the supplied WT/design/empirical sequence orientation
    relative to the reference. None permits inference only when unambiguous.
    """

    sample: str
    chrom: str
    donor_start: int  # 1-based inclusive
    donor_end: int  # 1-based inclusive
    donor_wt: str  # as provided (donor strand)
    donor_design: str
    donor_step2: str
    donor_strand: Optional[str] = None


@dataclass
class SampleOutcome:
    """Result of pre-calling one sample (mirrors one output line)."""

    sample: str
    editing_outcome: str
    designed_variants: List[str] = field(default_factory=list)
    synthetic_errors: List[str] = field(default_factory=list)
    spontaneous_variants: List[str] = field(default_factory=list)
    ambiguous_genotypes: List[dict] = field(default_factory=list)

    def to_row(self) -> str:
        return "\t".join(
            [
                self.sample,
                self.editing_outcome,
                ",".join(self.designed_variants),
                ",".join(self.synthetic_errors),
                ",".join(self.spontaneous_variants),
                json.dumps(self.ambiguous_genotypes, sort_keys=True, separators=(",", ":")),
            ]
        )


def orient_to_genome(spec: DonorSpec, chrom_seq: str) -> DonorSpec:
    """Orient all three donor sequences using explicit or unambiguous strand.

    The input specification is never mutated. A flipped result declares '+'
    because its sequences are now in reference orientation; the source declaration
    remains on the original specification and in the input samplesheet.

    Raises:
        ValueError: if neither orientation matches the reference window
            (the Perl ``exit 0``\\'s here; we raise so the driver can record it).
    """
    _validate_spec(spec)
    if spec.donor_end > len(chrom_seq):
        raise ValueError(f"donor window exceeds reference contig: {spec.sample}")
    window = chrom_seq[spec.donor_start - 1 : spec.donor_end].upper()
    plus_matches = spec.donor_wt == window
    minus_matches = reverse_complement(spec.donor_wt) == window
    if spec.donor_strand is not None:
        matches = plus_matches if spec.donor_strand == "+" else minus_matches
        if not matches:
            raise ValueError(f"declared donor_strand {spec.donor_strand!r} contradicts reference WT: {spec.sample}")
        flip = spec.donor_strand == "-"
    elif plus_matches and minus_matches:
        if any(seq != reverse_complement(seq) for seq in (spec.donor_design, spec.donor_step2)):
            raise ValueError(f"donor_strand is required for palindromic WT with asymmetric donor sequence: {spec.sample}")
        flip = False
    elif plus_matches or minus_matches:
        flip = minus_matches
    else:
        raise ValueError(f"provided donor sequence doesn't match reference genome: {spec.sample}")
    if not flip:
        return spec
    flipped = DonorSpec(
        sample=spec.sample,
        chrom=spec.chrom,
        donor_start=spec.donor_start,
        donor_end=spec.donor_end,
        donor_wt=reverse_complement(spec.donor_wt),
        donor_design=reverse_complement(spec.donor_design),
        donor_step2=reverse_complement(spec.donor_step2),
        donor_strand="+",
    )
    return flipped


@dataclass(frozen=True)
class Variant:
    """Structured genomic variant; the presentation ID is never parsed back."""

    chrom: str
    pos: int
    ref: str
    alt: str
    low_confidence: bool = False

    @property
    def identifier(self) -> str:
        return ("(LC)" if self.low_confidence else "") + f"{self.chrom}_{self.pos}_{self.ref}_to_{self.alt}"


def _vcf_rows(path: Path):
    """Read a single-sample (or eight-column site-only) VCF by named FORMAT keys."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        header_width = None
        for line_number, line in enumerate(fh, 1):
            if line.startswith("#CHROM\t"):
                if header_width is not None:
                    raise ValueError(f"{path}:{line_number}: repeated VCF header")
                header = line.rstrip("\r\n").split("\t")
                header_width = len(header)
                if (header_width not in (8, 10)
                        or header[:8] != ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
                        or (header_width == 10 and (header[8] != "FORMAT" or not header[9]))):
                    raise ValueError(f"{path}:{line_number}: invalid single-sample VCF header")
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\r\n").split("\t")
            if header_width is None or len(cols) != header_width:
                raise ValueError(f"{path}:{line_number}: missing header or mismatched VCF column count")
            try:
                pos = int(cols[1])
            except ValueError as exc:
                raise ValueError(f"{path}:{line_number}: invalid VCF POS {cols[1]!r}") from exc
            if pos < 1 or not cols[0] or not cols[3]:
                raise ValueError(f"{path}:{line_number}: invalid VCF coordinate or REF")
            fmt = {}
            alts = cols[4].upper().split(",")
            if any(not alt for alt in alts) or len(alts) != len(set(alts)) or ("." in alts and len(alts) > 1):
                raise ValueError(f"{path}:{line_number}: empty, duplicate, or mixed missing ALT")
            if len(cols) == 10:
                keys, values = cols[8].split(":"), cols[9].split(":")
                if any(not key for key in keys) or len(keys) != len(set(keys)) or len(values) > len(keys):
                    raise ValueError(f"{path}:{line_number}: invalid FORMAT fields")
                fmt = dict(zip(keys, values + ["."] * (len(keys) - len(values))))
                if fmt.get("AD", ".") != ".":
                    depths = fmt["AD"].split(",")
                    expected = 1 + (len(alts) if alts != ["."] else 0)
                    if len(depths) != expected or any(v != "." and (not v.isdigit() or int(v) < 0) for v in depths):
                        raise ValueError(f"{path}:{line_number}: AD must contain nonnegative depths for REF and every ALT")
            info = {}
            for item in cols[7].split(";") if cols[7] != "." else ():
                parts = item.split("=", 1)
                if not parts[0] or parts[0] in info:
                    raise ValueError(f"{path}:{line_number}: duplicate or blank INFO key")
                info[parts[0]] = parts[1] if len(parts) == 2 else None
            if "END" in info:
                try:
                    end = int(info["END"])
                except (ValueError, TypeError) as exc:
                    raise ValueError(f"{path}:{line_number}: END must be an integer") from exc
                if end < pos:
                    raise ValueError(f"{path}:{line_number}: END precedes POS")
            yield cols[0], pos, cols[3].upper(), alts, info, fmt
        if header_width is None:
            raise ValueError(f"{path}: no VCF #CHROM header")


def _selected_indices(alts: List[str], fmt: dict, *, site_only: bool = False) -> set[int]:
    """Validate genotype indexing, retaining unresolved selected alternatives."""
    if "GT" not in fmt:
        if site_only and len(alts) == 1:
            indices = {1} if alts != ["."] else set()
        elif all(a in (".", "<NON_REF>", "<*>") for a in alts):
            return set()
        else:
            raise ValueError("concrete ALT requires GT (or a biallelic site-only VCF)")
    else:
        try:
            indices = {int(a) for a in fmt["GT"].replace("|", "/").split("/") if a != "."}
        except ValueError as exc:
            raise ValueError(f"invalid GT {fmt['GT']!r}") from exc
        if any(i < 0 or i > len(alts) for i in indices):
            raise ValueError(f"GT allele outside ALT list: {fmt['GT']!r}")
        indices.discard(0)
    if any(alts[i - 1] == "." for i in indices):
        raise ValueError("GT selects missing ALT")
    return indices


def _selected_alts(alts: List[str], fmt: dict, *, site_only: bool = False) -> List[str]:
    """Return concrete clonal alleles; raw unresolved evidence is recorded separately."""
    indices = _selected_indices(alts, fmt, site_only=site_only)
    if _ambiguous_genotype(fmt):
        return []
    return [alts[i - 1] for i in sorted(indices) if alts[i - 1] and not set(alts[i - 1]) - set("ACGT")]


def _ambiguous_genotype(fmt: dict) -> bool:
    if "GT" not in fmt:
        return False
    values = fmt["GT"].replace("|", "/").split("/")
    return "." in values or len(set(values)) > 1


def _iter_vcf_records(path: Path):
    """Yield genotype-selected sequence variants; FILTER remains legacy-compatible."""
    for chrom, pos, ref, alts, _info, fmt in _vcf_rows(path):
        for alt in _selected_alts(alts, fmt, site_only=not fmt):
            yield chrom, pos, ref, alt


def _reference_depth(fmt: dict, *, reference_block: bool) -> int:
    """Conservative reference depth; missing values never imply coverage.

    For hom-reference blocks MIN_DP is a lower bound over the whole END interval.
    AD/DP apply only to an explicit position; a block lacking MIN_DP cannot prove
    every interior base covered.
    """
    if reference_block:
        value = fmt.get("MIN_DP", ".")
    elif "AD" in fmt and fmt["AD"] != ".":
        value = fmt["AD"].split(",")[0]
    elif fmt.get("GT", ".").replace("|", "/") in ("0", "0/0"):
        value = fmt.get("DP", ".")
    else:
        value = "."
    if value in ("", "."):
        return 0
    try:
        depth = int(value)
    except ValueError as exc:
        raise ValueError(f"invalid reference depth {value!r}") from exc
    if depth < 0:
        raise ValueError(f"negative reference depth {depth}")
    return depth


def precall_sample(
    spec: DonorSpec,
    vcf_path: Path,
    gvcf_path: Path,
    chrom_seq: str,
    ne_threshold: int = 1,
    quantification_window_size: int = 10,
    quantification_window_cov_perc: float = 1.0,
) -> SampleOutcome:
    """Pre-call the on-target editing outcome for a single sample.

    Classify a donor using genotype-selected alleles and conservative reference coverage.
    """
    if ne_threshold < 1 or quantification_window_size < 1:
        raise ValueError("coverage threshold and window size must be positive")
    if not 0 < quantification_window_cov_perc <= 1:
        raise ValueError("coverage fraction must be in (0, 1]")
    spec = orient_to_genome(spec, chrom_seq)
    donor_wt, donor_design, donor_step2 = spec.donor_wt, spec.donor_design, spec.donor_step2
    dstart, dend = spec.donor_start, spec.donor_end
    selected: Dict[Tuple[str, int, str, str], Tuple[Variant, str]] = {}
    called_spans: List[Tuple[int, int]] = []
    ambiguous: List[dict] = []
    unresolved = False

    def note_genotype(chrom, pos, ref, alts, fmt, *, site_only=False, span_end=None):
        nonlocal unresolved
        span_end = pos + len(ref) - 1 if span_end is None else span_end
        if chrom != spec.chrom or pos > dend or span_end < dstart:
            return
        if _ambiguous_genotype(fmt):
            evidence = {"chrom": chrom, "pos": pos, "GT": fmt["GT"]}
            if evidence not in ambiguous:
                ambiguous.append(evidence)
        indices = _selected_indices(alts, fmt, site_only=site_only)
        if any(not alts[i - 1] or set(alts[i - 1]) - set("ACGT") for i in indices):
            unresolved = True
            evidence = {"chrom": chrom, "pos": pos, "GT": fmt.get("GT", "."),
                        "ALT": ",".join(alts), "reason": "unresolved_allele"}
            if evidence not in ambiguous:
                ambiguous.append(evidence)

    def add_variant(chrom, pos, ref, alt, *, low_confidence=False):
        pos0 = pos - dstart
        if chrom != spec.chrom:
            return
        if pos0 < 0 or pos0 + len(ref) > len(donor_wt):
            if pos <= dend and pos + len(ref) - 1 >= dstart:
                raise ValueError(f"variant overlaps donor boundary and cannot be represented: {chrom}:{pos}")
            return
        key = (chrom, pos, ref, alt)
        if key in selected:
            return  # exact duplicate is one observation, not extra support
        end = pos + len(ref) - 1
        if low_confidence and any(pos <= hi and end >= lo for lo, hi in called_spans):
            return
        for old, _bucket in selected.values():
            if pos <= old.pos + len(old.ref) - 1 and end >= old.pos:
                raise ValueError(f"overlapping/conflicting donor alleles: {old.identifier}, {key}")
        donor_tmp = apply_variant(donor_wt, pos0, ref, alt)
        bucket, _e2, _e1 = bucket_variant(donor_tmp, donor_wt, donor_design, donor_step2)
        if not low_confidence:
            called_spans.append((pos, end))
        if low_confidence and bucket == SPONTANEOUS:
            return  # preserve the documented no-spontaneous-LC convention
        selected[key] = (Variant(chrom, pos, ref, alt, low_confidence), bucket)

    for chrom, pos, ref, alts, info, fmt in _vcf_rows(vcf_path):
        span_end = max(pos + len(ref) - 1, int(info.get("END", pos)))
        if chrom != spec.chrom or pos > dend or span_end < dstart:
            continue
        chosen = _selected_alts(alts, fmt, site_only=not fmt)
        note_genotype(chrom, pos, ref, alts, fmt, site_only=not fmt, span_end=span_end)
        for alt in chosen:
            add_variant(chrom, pos, ref, alt)

    # The actual clipped donor window supplies the denominator, including short donors.
    width = min(quantification_window_size, len(donor_wt))
    qstart = math.ceil((dstart + dend - width) / 2)
    qend = qstart + width - 1
    covered: set[int] = set()
    # Coverage at the donor center cannot establish absence of an edit elsewhere.
    # Same-length designs have explicit positional differences. A length change
    # can align in several ways (especially in repeats), so require the entire
    # donor interval rather than invent a single diagnostic coordinate.
    diagnostic = ({dstart + i for i, (wt, design) in enumerate(zip(donor_wt, donor_design)) if wt != design}
                  if len(donor_wt) == len(donor_design) else set(range(dstart, dend + 1)))
    seen_gvcf = []
    for chrom, pos, ref, alts, info, fmt in _vcf_rows(gvcf_path):
        if chrom != spec.chrom:
            continue
        if pos + len(ref) - 1 > len(chrom_seq) or chrom_seq[pos - 1:pos - 1 + len(ref)].upper() != ref:
            raise ValueError(f"gVCF REF does not match reference at {chrom}:{pos}")
        end = int(info.get("END", pos))
        if end < pos or end > len(chrom_seq):
            raise ValueError(f"invalid gVCF END at {chrom}:{pos}: {end}")
        symbolic = all(a in (".", "<NON_REF>", "<*>") for a in alts)
        hom_ref = fmt.get("GT", ".").replace("|", "/") in ("0", "0/0")
        is_block = end > pos and symbolic and hom_ref
        if end > pos and not symbolic and all(not set(a) - set("ACGTN") for a in alts):
            raise ValueError(f"END on concrete sequence allele is not a reference block: {chrom}:{pos}")
        span_end = max(end, pos + len(ref) - 1)
        if pos <= dend and span_end >= dstart:
            identity = (pos, span_end, ref, tuple(alts), tuple(sorted(info.items())), tuple(sorted(fmt.items())))
            if identity in seen_gvcf:
                continue
            # planning/328: with -ERC BP_RESOLUTION, GATK still emits a single-base <NON_REF> record at every base
            # a called deletion's REF spans (ATC->A at 995605, then T at 995606 and C at 995607). Those records
            # restate bases the deletion record already describes; they are not a second, conflicting call. A
            # dense MNV design is exactly what GATK writes as such an insertion/deletion pair. Skipped entirely:
            # they add no variant (symbolic ALT) and no reference coverage that the deletion has not settled.
            if (symbolic and len(ref) == 1 and end == pos
                    and any(previous[0] < pos <= previous[1] and len(previous[2]) > 1
                            and any(a not in (".", "<NON_REF>", "<*>") for a in previous[3])
                            for previous in seen_gvcf)):
                continue
            if any(pos <= previous[1] and span_end >= previous[0] for previous in seen_gvcf):
                raise ValueError(f"overlapping non-identical gVCF records in donor interval: {chrom}:{pos}")
            seen_gvcf.append(identity)
        depth = _reference_depth(fmt, reference_block=is_block)
        if depth >= ne_threshold:
            coverage_end = end if is_block else pos
            covered.update(range(max(pos, dstart), min(coverage_end, dend) + 1))
        for alt in _selected_alts(alts, fmt):
            add_variant(chrom, pos, ref, alt, low_confidence=True)
        note_genotype(chrom, pos, ref, alts, fmt, span_end=span_end)

    designed = [v for v, bucket in selected.values() if bucket == DESIGNED]
    synthetic = [v for v, bucket in selected.values() if bucket == SYNTHETIC_ERROR]
    spontaneous = [v for v, bucket in selected.values() if bucket == SPONTANEOUS]

    def reconstruct(variants):
        result = donor_wt
        # Descending original-reference coordinates preserve offsets across indels.
        for variant in sorted(variants, key=lambda v: v.pos, reverse=True):
            result = apply_variant(result, variant.pos - dstart, variant.ref, variant.alt)
        return result

    combined = reconstruct(designed)
    combined_all = reconstruct(designed + spontaneous)
    central_covered = sum(qstart <= pos <= qend for pos in covered)
    proves_unedited = central_covered >= width * quantification_window_cov_perc and diagnostic <= covered
    if donor_design == donor_wt:
        flag = "designed_variant_unedited" if proves_unedited else "unknown"
    elif combined == donor_design:
        flag = "designed_variant_installed"
    elif combined_all == donor_design:
        flag = "designed_variant_installed"
        designed.extend(spontaneous)
        spontaneous = []
    elif levenshtein(combined, donor_design) < levenshtein(donor_wt, donor_design):
        flag = "designed_variant_partially_installed"
    elif proves_unedited:
        flag = "designed_variant_unedited"
    else:
        flag = "unknown"
    if synthetic:
        flag += ";synthetic_error_installed"
    if spontaneous:
        flag += ";spontaneous_mutations"
    if ambiguous:
        reasons = []
        if any(e.get("reason") != "unresolved_allele" for e in ambiguous):
            reasons.append("ambiguous_genotype")
        if unresolved:
            reasons.append("unresolved_allele")
        flag = "unknown;" + ";".join(reasons)
    return SampleOutcome(
        spec.sample, flag,
        sorted(v.identifier for v in designed),
        sorted(v.identifier for v in synthetic),
        sorted(v.identifier for v in spontaneous),
        ambiguous,
    )


# --------------------------------------------------------------------------- #
# Samplesheet + driver (reproduces the Perl CLI contract)                      #
# --------------------------------------------------------------------------- #

# Canonical samplesheet column names (== Shengdi's positional columns 0,3-8)
_COL = {
    "sample": "sample",
    "chrom": "chrom",
    "donor_start": "donor_start_coord",
    "donor_end": "donor_end_coord",
    "donor_wt": "wt_seq",
    "donor_design": "design_donor",
    "donor_step2": "synthetic_donor",
}


def _validate_spec(spec: DonorSpec) -> None:
    if spec.donor_strand not in (None, "+", "-"):
        raise ValueError(f"donor_strand must be '+' or '-' when supplied: {spec.sample}")
    if not spec.sample or not spec.chrom or any(c in spec.sample for c in "/\\\t\n\r"):
        raise ValueError("sample and contig names must be nonempty; sample must be a safe filename")
    if spec.sample in (".", "..") or spec.donor_start < 1 or spec.donor_end < spec.donor_start:
        raise ValueError(f"invalid donor coordinates or sample: {spec.sample}")
    if len(spec.donor_wt) != spec.donor_end - spec.donor_start + 1:
        raise ValueError(f"wt_seq length differs from 1-based inclusive donor window: {spec.sample}")
    for label, value in (("wt_seq", spec.donor_wt), ("design_donor", spec.donor_design),
                         ("synthetic_donor", spec.donor_step2)):
        if not value or set(value) - set("ACGTN"):
            raise ValueError(f"invalid {label} DNA for {spec.sample}")


def read_samplesheet(path: Path) -> List[DonorSpec]:
    """Read and validate a unique-sample, 1-based inclusive donor samplesheet."""
    specs: List[DonorSpec] = []
    seen = set()
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        header = reader.fieldnames or []
        if len(header) != len(set(header)) or set(_COL.values()) - set(header):
            raise ValueError(f"{path}: missing or duplicate samplesheet columns")
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f"{path}:{reader.line_num}: malformed samplesheet row")
            fields = {key: row[col].strip() for key, col in _COL.items()}
            spec = DonorSpec(fields["sample"], fields["chrom"],
                             int(fields["donor_start"]), int(fields["donor_end"]),
                             fields["donor_wt"].upper(), fields["donor_design"].upper(),
                             fields["donor_step2"].upper(),
                             row.get("donor_strand", "").strip() or None)
            _validate_spec(spec)
            if spec.sample in seen:
                raise ValueError(f"{path}: duplicate sample {spec.sample!r}")
            seen.add(spec.sample)
            specs.append(spec)
    if not specs:
        raise ValueError(f"{path}: samplesheet contains no samples")
    return specs


def precall(
    samplesheet: Path,
    workdir: Path,
    genome: Path,
    out: Path,
    ne_threshold: int = 1,
    quantification_window_size: int = 10,
    quantification_window_cov_perc: float = 1.0,
    samples: Optional[set] = None,
    *,
    loaded_reference: Optional[LoadedReference] = None,
) -> List[SampleOutcome]:
    """Run pre-calling over a samplesheet; write the Shengdi-format summary.

    Reads ``{workdir}/vcf/{sample}.fil.vcf`` and ``{workdir}/gvcf/{sample}.g.vcf``.
    If ``samples`` is given, only those sample names are processed (for parity
    subsets).
    """
    refs = loaded_reference.require_current(genome) if loaded_reference is not None else load_fasta(str(genome))
    specs = read_samplesheet(samplesheet)
    if samples is not None:
        absent = set(samples) - {spec.sample for spec in specs}
        if absent:
            raise ValueError(f"requested samples absent from samplesheet: {sorted(absent)}")
        specs = [spec for spec in specs if spec.sample in samples]
    if not specs:
        raise ValueError("no samples selected for precall")
    missing = [str(workdir / subdir / f"{spec.sample}{suffix}")
               for spec in specs for subdir, suffix in (("vcf", ".fil.vcf"), ("gvcf", ".g.vcf"))
               if not (workdir / subdir / f"{spec.sample}{suffix}").is_file()]
    if missing:
        raise FileNotFoundError("precall missing required input(s): " + ", ".join(missing))
    for spec in specs:
        if spec.chrom not in refs:
            raise ValueError(f"sample {spec.sample}: donor contig {spec.chrom!r} absent from reference")
        orient_to_genome(spec, refs[spec.chrom])
    results: List[SampleOutcome] = []
    out.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=out.name + ".", suffix=".tmp", dir=out.parent)
    try:
        with os.fdopen(fd, "w") as oh:
            oh.write(
                "#sample_name\tediting_outcome\tdesigned_variants\t"
                "synthetic_errors\tspontaneous_variants\tambiguous_genotypes\n"
            )
            for spec in specs:
                outcome = precall_sample(
                    spec, workdir / "vcf" / f"{spec.sample}.fil.vcf",
                    workdir / "gvcf" / f"{spec.sample}.g.vcf", refs[spec.chrom],
                    ne_threshold=ne_threshold,
                    quantification_window_size=quantification_window_size,
                    quantification_window_cov_perc=quantification_window_cov_perc,
                )
                results.append(outcome)
                oh.write(outcome.to_row() + "\n")
        os.replace(temporary, out)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)
    return results


def _build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Pre-call on-target editing outcomes from WGS VCF/gVCF "
        "(Python port of Shengdi Li's precalling_target_outcome.pl)."
    )
    p.add_argument("--in", dest="samplesheet", required=True, type=Path)
    p.add_argument("--workdir", required=True, type=Path)
    p.add_argument("--genome", required=True, type=Path)
    p.add_argument("--out", required=True, type=Path)
    p.add_argument("--NE-threshold", dest="ne_threshold", type=int, default=1)
    p.add_argument("--quantification-window-size", type=int, default=10)
    p.add_argument("--quantification-window-cov-perc", type=float, default=1.0)
    return p


def main(argv: Optional[List[str]] = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    precall(
        samplesheet=args.samplesheet,
        workdir=args.workdir,
        genome=args.genome,
        out=args.out,
        ne_threshold=args.ne_threshold,
        quantification_window_size=args.quantification_window_size,
        quantification_window_cov_perc=args.quantification_window_cov_perc,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
