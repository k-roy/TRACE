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

The port reproduces the Perl semantics exactly, including: no FILTER-column
filtering of the ``.fil.vcf`` (lowQual in-window variants are still bucketed);
strict ``<`` comparisons; gVCF AD taken as the 2nd FORMAT field with ref depth =
its first value; the donor-region "combine" reconstruction used to separate
complete from partial installation.
"""

from __future__ import annotations

import argparse
import math
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
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    refs[name] = "".join(chunks).upper()
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        refs[name] = "".join(chunks).upper()
    return refs


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
    return donor_wt[:pos0] + alt + donor_wt[pos0 + len(ref) :]


# --------------------------------------------------------------------------- #
# Per-sample pre-calling                                                       #
# --------------------------------------------------------------------------- #


@dataclass
class DonorSpec:
    """A sample's donor-window specification (one samplesheet row)."""

    sample: str
    chrom: str
    donor_start: int  # 1-based inclusive
    donor_end: int  # 1-based inclusive
    donor_wt: str  # as provided (donor strand)
    donor_design: str
    donor_step2: str


@dataclass
class SampleOutcome:
    """Result of pre-calling one sample (mirrors one output line)."""

    sample: str
    editing_outcome: str
    designed_variants: List[str] = field(default_factory=list)
    synthetic_errors: List[str] = field(default_factory=list)
    spontaneous_variants: List[str] = field(default_factory=list)

    def to_row(self) -> str:
        return "\t".join(
            [
                self.sample,
                self.editing_outcome,
                ",".join(self.designed_variants),
                ",".join(self.synthetic_errors),
                ",".join(self.spontaneous_variants),
            ]
        )


def orient_to_genome(spec: DonorSpec, chrom_seq: str) -> DonorSpec:
    """Flip donors to genome orientation if the wt donor is on the minus strand.

    Raises:
        ValueError: if neither orientation matches the reference window
            (the Perl ``exit 0``\\'s here; we raise so the driver can record it).
    """
    window = chrom_seq[spec.donor_start - 1 : spec.donor_end]
    if spec.donor_wt == window:
        return spec
    flipped = DonorSpec(
        sample=spec.sample,
        chrom=spec.chrom,
        donor_start=spec.donor_start,
        donor_end=spec.donor_end,
        donor_wt=reverse_complement(spec.donor_wt),
        donor_design=reverse_complement(spec.donor_design),
        donor_step2=reverse_complement(spec.donor_step2),
    )
    if flipped.donor_wt == window:
        return flipped
    raise ValueError(f"provided donor sequence doesn't match reference genome: {spec.sample}")


def _iter_vcf_records(path: Path):
    """Yield ``(chrom, pos, ref, alt_field)`` for non-header VCF lines."""
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            # chrom, pos, id, ref, alt  (FILTER intentionally ignored, per Perl)
            yield cols[0], int(cols[1]), cols[3], cols[4]


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

    Faithful port of ``precalling_target_outcome.pl`` per-sample logic.
    """
    spec = orient_to_genome(spec, chrom_seq)
    donor_wt = spec.donor_wt
    donor_design = spec.donor_design
    donor_step2 = spec.donor_step2
    chrom = spec.chrom
    dstart, dend = spec.donor_start, spec.donor_end
    offset = dstart - 1
    wt_len = len(donor_wt)

    designed: Dict[str, str] = {}
    synthetic: Dict[str, str] = {}
    spontaneous: Dict[str, str] = {}
    called_positions: set = set()  # donor-local 0-based positions called in VCF

    # ---- VCF pass (all non-header records; NO FILTER filtering, per Perl) ----
    for v_chr, v_pos, v_ref, v_alt in _iter_vcf_records(vcf_path):
        pos0 = v_pos - offset - 1
        if v_chr != chrom or pos0 < 0 or (pos0 + len(v_ref)) >= wt_len:
            continue
        donor_tmp = apply_variant(donor_wt, pos0, v_ref, v_alt)
        for k in range(len(v_ref)):
            called_positions.add(pos0 + k)
        bucket, e2, e1 = bucket_variant(donor_tmp, donor_wt, donor_design, donor_step2)
        v_id = f"{v_chr}_{v_pos}_{v_ref}_to_{v_alt}"
        evidence = f"{e2}->{e1}"
        if bucket == DESIGNED:
            designed[v_id] = evidence
        elif bucket == SYNTHETIC_ERROR:
            synthetic[v_id] = evidence
        else:
            spontaneous[v_id] = evidence

    # ---- gVCF pass: coverage check + less-confident (LC) promotions ----------
    flag_covered = 0
    qwin_start = math.ceil(dend - (dend - dstart) / 2 - quantification_window_size / 2)
    qwin_end = qwin_start + quantification_window_size - 1
    with open(gvcf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            v_chr, v_pos = cols[0], int(cols[1])
            v_ref, v_alts = cols[3], cols[4]
            v_alt = v_alts.split(",")
            # AD = 2nd FORMAT value; ref depth = first AD entry
            sample_field = cols[9]
            info = sample_field.split(":")
            ref_depth = int(info[1].split(",")[0]) if len(info) > 1 else 0

            if v_chr == chrom and ref_depth >= ne_threshold and qwin_start <= v_pos <= qwin_end:
                flag_covered += 1

            if v_chr != chrom:
                continue
            if not (offset + 1 <= v_pos < offset + 1 + wt_len - len(v_ref)):
                continue
            if v_alt[0] == "<NON_REF>":
                continue

            pos0 = v_pos - offset - 1
            donor_tmp = apply_variant(donor_wt, pos0, v_ref, v_alt[0])
            bucket, e2, e1 = bucket_variant(donor_tmp, donor_wt, donor_design, donor_step2)
            already_called = any((pos0 + k) in called_positions for k in range(len(v_ref)))
            if already_called:
                continue
            v_id = f"{v_chr}_{v_pos}_{v_ref}_to_{v_alt[0]}"
            evidence = f"{e2}->{e1}"
            if bucket == DESIGNED:
                designed[f"(LC){v_id}"] = evidence
            elif bucket == SYNTHETIC_ERROR:
                synthetic[f"(LC){v_id}"] = evidence
            # NB: spontaneous LC calls are intentionally NOT recorded (per Perl)

    # ---- combine designed (and designed+spontaneous) to judge installation ---
    combine = list(donor_wt)
    combine_all = list(donor_wt)
    for v_id in designed:
        sp = v_id.split("_")
        start_pos = int(sp[1]) - offset - 1
        ref_len = len(sp[2])
        for k in range(ref_len):
            combine[start_pos + k] = ""
            combine_all[start_pos + k] = ""
        combine[start_pos] = sp[4]
        combine_all[start_pos] = sp[4]
    for v_id in spontaneous:
        sp = v_id.split("_")
        start_pos = int(sp[1]) - offset - 1
        ref_len = len(sp[2])
        for k in range(ref_len):
            combine_all[start_pos + k] = ""
        combine_all[start_pos] = sp[4]
    donor_region_combine = "".join(combine)
    donor_region_combine_all = "".join(combine_all)

    if levenshtein(donor_region_combine, donor_design) == 0:
        flag = "designed_variant_installed"
    elif levenshtein(donor_region_combine_all, donor_design) == 0:
        flag = "designed_variant_installed"
        # designed + spontaneous reconstruct the design: promote spontaneous
        for k, v in spontaneous.items():
            designed[k] = v
        spontaneous = {}
    elif levenshtein(donor_region_combine, donor_design) < levenshtein(donor_wt, donor_design):
        flag = "designed_variant_partially_installed"
    elif flag_covered >= quantification_window_size * quantification_window_cov_perc:
        flag = "designed_variant_unedited"
    else:
        flag = "unknown"

    editing_outcome = flag
    if synthetic:
        editing_outcome += ";synthetic_error_installed"
    if spontaneous:
        editing_outcome += ";spontaneous_mutations"

    return SampleOutcome(
        sample=spec.sample,
        editing_outcome=editing_outcome,
        designed_variants=sorted(designed),
        synthetic_errors=sorted(synthetic),
        spontaneous_variants=sorted(spontaneous),
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


def read_samplesheet(path: Path) -> List[DonorSpec]:
    """Read a tab-separated samplesheet (with header) into :class:`DonorSpec`\\ s."""
    specs: List[DonorSpec] = []
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {name: header.index(col) for name, col in _COL.items()}
        for line in fh:
            if not line.strip():
                continue
            c = line.rstrip("\n").split("\t")
            specs.append(
                DonorSpec(
                    sample=c[idx["sample"]],
                    chrom=c[idx["chrom"]],
                    donor_start=int(c[idx["donor_start"]]),
                    donor_end=int(c[idx["donor_end"]]),
                    donor_wt=c[idx["donor_wt"]].upper(),
                    donor_design=c[idx["donor_design"]].upper(),
                    donor_step2=c[idx["donor_step2"]].upper(),
                )
            )
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
) -> List[SampleOutcome]:
    """Run pre-calling over a samplesheet; write the Shengdi-format summary.

    Reads ``{workdir}/vcf/{sample}.fil.vcf`` and ``{workdir}/gvcf/{sample}.g.vcf``.
    If ``samples`` is given, only those sample names are processed (for parity
    subsets).
    """
    refs = load_fasta(str(genome))
    specs = read_samplesheet(samplesheet)
    results: List[SampleOutcome] = []
    with open(out, "w") as oh:
        oh.write(
            "#sample_name\tediting_outcome\tdesigned_variants\t"
            "synthetic_errors\tspontaneous_variants\n"
        )
        for spec in specs:
            if samples is not None and spec.sample not in samples:
                continue
            vcf_path = workdir / "vcf" / f"{spec.sample}.fil.vcf"
            gvcf_path = workdir / "gvcf" / f"{spec.sample}.g.vcf"
            if not vcf_path.exists() or not gvcf_path.exists():
                continue
            outcome = precall_sample(
                spec,
                vcf_path,
                gvcf_path,
                refs[spec.chrom],
                ne_threshold=ne_threshold,
                quantification_window_size=quantification_window_size,
                quantification_window_cov_perc=quantification_window_cov_perc,
            )
            results.append(outcome)
            oh.write(outcome.to_row() + "\n")
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
