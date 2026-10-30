"""
trace_crispr.wgs.target_calls
=============================

Haplotype-level editing calls for one site, with the target codon judged apart from its synonymous helpers.

Rules ruled by Kevin R. Roy, 2026-10-30 (SGTC planning/328l, items 2.1-2.3), on read pairs (fragments) whose
designed bases were read by local realignment:

* **2.1, read-pair floors from contamination.** ``c`` is the chance that a read pair comes from another
  genome (measured at the bc1 locus: 0.6% other samples + about 5% other cells in the colony, so 0.054 by
  default). One read pair carrying the design calls it present; no edit needs >= 2 wild-type read pairs and
  none carrying the design; a mixture is called only when the minority wild-type count is improbable under
  ``c`` (binomial upper tail < ``alpha``, 0.01), otherwise the majority call stands.
* **2.2.1** a call rests on the designed bases that were read; the rest are flagged, not required.
* **2.2.2** "I would call 2.2.2 an incomplete edit when a designed base has a wt position. If those are in the
  syn recoded portion, then it's fine to have a separate category that defines correct aa seq. In general,
  for amino acid based edits, we should score editing outcomes based on the resulting translation frame of
  each gene." And: "Syn designs can be treated liked missense for the syn helper mutations, which are distint
  from the target codon." And: "separate category, helpers only."

Two calls come out:

``call_nt`` (nucleotide level): IE if any read pair carries the design at some designed bases and wild type at
others; else CE / MIX / NE / UN by the 2.1 floors.

``call_target`` (target level): each designed base is placed in its codon of the CDS (bases outside it are
protein-neutral helpers). A codon is read only when all its designed bases were read in that read pair. The
target is judged on its amino acid for a missense design and on its codon for a synonymous one; a helper codon
reading a third amino acid (a half-installed helper making a missense) counts as "other".

* ``CE``            target as designed, every designed base read as design
* ``CE_target``     target as designed, some helper read wild type
* ``helpers_only``  target wild type in read pairs whose helpers carry the design
* ``IE``            another amino acid or codon at the target or a helper
* ``MIX``           the minority genotype is more than ``c`` explains, or target-design and helpers-only
                    molecules coexist
* ``NE``            >= 2 read pairs read the target as wild type, with wild-type helpers or none read
* ``CE_unconfirmed_target``  design read only at helper bases; no read pair reached the target
* ``UN``            otherwise

A base that is neither the reference nor the designed one (``O``) is a third allele or an error, not a partial
edit; such read pairs make IE only when wild-type read pairs do not outnumber them and none carries the design.

Author: Kevin R. Roy
"""

from __future__ import annotations

import functools
from dataclasses import dataclass
from math import comb
from typing import Dict, FrozenSet, Iterable, List, Optional, Sequence, Tuple

DEFAULT_CONTAMINATION = 0.054
DEFAULT_MINORITY_ALPHA = 0.01

CE, CE_TARGET, HELPERS_ONLY, IE, MIX, NE, UN = "CE", "CE_target", "helpers_only", "IE", "MIX", "NE", "UN"
CE_UNCONFIRMED_TARGET = "CE_unconfirmed_target"

_BASES = "TCAG"
_AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
GENETIC_CODE = {a + b + c: _AAS[16 * i + 4 * j + k] for i, a in enumerate(_BASES) for j, b in enumerate(_BASES)
                for k, c in enumerate(_BASES)}
_COMP = str.maketrans("ACGT", "TGCA")


def _rc(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


@dataclass(frozen=True)
class CDS:
    """A single-exon coding sequence in the SAMPLE reference's coordinates (1-based, inclusive)."""

    chrom: str
    start: int
    end: int
    strand: str

    def __post_init__(self):
        if self.strand not in "+-" or len(self.strand) != 1 or self.start < 1 or self.end < self.start:
            raise ValueError(f"invalid CDS {self}")
        if (self.end - self.start + 1) % 3:
            raise ValueError(f"CDS length is not a multiple of 3: {self}")

    def codon_of(self, pos: int) -> Optional[Tuple[int, int]]:
        """(codon index from the start codon, offset within the codon in coding orientation), or None."""
        if not self.start <= pos <= self.end:
            return None
        d = pos - self.start if self.strand == "+" else self.end - pos
        return d // 3, d % 3

    def codon(self, fasta, k: int) -> str:
        if self.strand == "+":
            return fasta.fetch(self.chrom, self.start - 1 + 3 * k, self.start - 1 + 3 * k + 3).upper()
        return _rc(fasta.fetch(self.chrom, self.end - 3 * k - 3, self.end - 3 * k).upper())


@dataclass(frozen=True)
class TargetSpec:
    """Where a design's target codon sits: its CDS and the codon's genomic span (1-based, inclusive)."""

    cds: CDS
    codon_start: int
    codon_end: int


@dataclass(frozen=True)
class TargetModel:
    """A site's designed bases mapped onto codons, with the target separated from the helpers."""

    n: int
    codons: Dict[int, Tuple[Tuple[int, int, str], ...]]   # codon -> (allele index, offset, coding-strand alt)
    target: int
    helpers: FrozenSet[int]                               # allele indices outside the target codon
    wt_codon: Dict[int, str]
    wt_aa: Dict[int, str]
    design_aa: Dict[int, str]
    target_wt_codon: str
    target_design_codon: str
    synonymous_target: bool


@functools.lru_cache(maxsize=64)
def _verified_cds(fasta_path: str, cds: CDS) -> bool:
    import pysam
    with pysam.FastaFile(fasta_path) as fa:
        seq = fa.fetch(cds.chrom, cds.start - 1, cds.end).upper()
    seq = seq if cds.strand == "+" else _rc(seq)
    prot = "".join(GENETIC_CODE.get(seq[i:i + 3], "X") for i in range(0, len(seq), 3))
    if not (prot[0] == "M" and prot[-1] == "*" and "*" not in prot[:-1] and "X" not in prot):
        raise ValueError(f"{cds} does not translate M...* on {fasta_path}: {prot[:5]}...{prot[-5:]}")
    return True


def build_target_model(fasta, alleles: Sequence, spec: TargetSpec) -> TargetModel:
    """Map single-base designed substitutions onto `spec.cds` codons. Refuses anything it cannot verify:
    a CDS that does not translate M...*, a target codon off-frame or outside the CDS, a REF mismatch."""
    cds = spec.cds
    _verified_cds(fasta.filename.decode() if isinstance(fasta.filename, bytes) else str(fasta.filename), cds)
    t0, t1 = cds.codon_of(spec.codon_start), cds.codon_of(spec.codon_end)
    if t0 is None or t1 is None or t0[0] != t1[0] or spec.codon_end - spec.codon_start != 2:
        raise ValueError(f"target codon {spec.codon_start}-{spec.codon_end} is not one codon of {cds}")
    tk = t0[0]
    codons: Dict[int, List[Tuple[int, int, str]]] = {}
    for i, a in enumerate(alleles):
        if len(a.ref) != 1 or len(a.alt) != 1:
            raise ValueError("target-level calling supports single-base designed substitutions only")
        if fasta.fetch(cds.chrom, a.pos - 1, a.pos).upper() != a.ref:
            raise ValueError(f"{cds.chrom}:{a.pos}: REF {a.ref} does not match the reference")
        co = cds.codon_of(a.pos)
        if co is not None:
            codons.setdefault(co[0], []).append((i, co[1], a.alt if cds.strand == "+" else a.alt.translate(_COMP)))
    if tk not in codons:
        raise ValueError("no designed base falls in the target codon")
    wt_codon = {k: cds.codon(fasta, k) for k in codons}
    design_aa = {}
    for k, items in codons.items():
        c = list(wt_codon[k])
        for _i, off, b in items:
            c[off] = b
        design_aa[k] = GENETIC_CODE["".join(c)]
    target_design = list(wt_codon[tk])
    for _i, off, b in codons[tk]:
        target_design[off] = b
    target_design = "".join(target_design)
    helpers = frozenset(i for i in range(len(alleles)) if i not in {j for j, _, _ in codons[tk]})
    return TargetModel(len(alleles), {k: tuple(v) for k, v in codons.items()}, tk, helpers, wt_codon,
                       {k: GENETIC_CODE[c] for k, c in wt_codon.items()}, design_aa, wt_codon[tk], target_design,
                       GENETIC_CODE[wt_codon[tk]] == GENETIC_CODE[target_design])


def binom_tail(k: int, n: int, p: float) -> float:
    """P(X >= k) for X ~ Binomial(n, p)."""
    return sum(comb(n, i) * p ** i * (1 - p) ** (n - i) for i in range(k, n + 1))


def significant_minority(n_major: int, n_minor: int, c: float, alpha: float) -> bool:
    return n_minor > 0 and binom_tail(n_minor, n_major + n_minor, c) < alpha


def fragment_haplotypes(observations: Iterable[Sequence[Optional[str]]]) -> List[str]:
    """Realigner observations per fragment ("alt", "ref", another base, "conflict", None) -> 'D', 'W', 'O', '.'."""
    out = []
    for obs in observations:
        out.append("".join("D" if o == "alt" else "W" if o == "ref" else "." if o in (None, "conflict") else "O"
                           for o in obs))
    return out


@dataclass(frozen=True)
class TargetCall:
    call_nt: str
    call_target: str
    n_design: int
    n_wt: int
    n_partial: int
    n_other: int
    target_design: int
    target_wt: int
    target_helpers_only: int
    target_other: int
    target_unexpected_base: int
    helpers_design_target_unread: int
    incomplete_coverage: bool
    synonymous_target: bool


def call_haplotypes(frags: Sequence[str], model: TargetModel, *, c: float = DEFAULT_CONTAMINATION,
                    alpha: float = DEFAULT_MINORITY_ALPHA) -> TargetCall:
    """Both calls for one site from its fragment haplotypes (strings of D/W/O/. in allele order)."""
    if not 0 < c < 1 or not 0 < alpha < 1:
        raise ValueError("contamination rate and alpha must lie in (0, 1)")
    nD = nW = nP = nO = 0
    pD = pW = pX = pXo = pH = uH = 0
    helpers_wt_in_design = False
    read_d, read_w = set(), set()
    tk = model.target
    tgt = model.codons[tk]
    for cs in frags:
        if len(cs) != model.n:
            raise ValueError("fragment haplotype length does not match the site's alleles")
        hasd, hasw, haso = "D" in cs, "W" in cs, "O" in cs
        if haso:
            nO += 1
        elif hasd and hasw:
            nP += 1
        elif hasd:
            nD += 1
        elif hasw:
            nW += 1
        read_d |= {i for i, x in enumerate(cs) if x == "D"}
        read_w |= {i for i, x in enumerate(cs) if x == "W"}
        states = {}
        for k, items in model.codons.items():
            if any(cs[i] == "." for i, _, _ in items):
                continue
            if any(cs[i] == "O" for i, _, _ in items):
                states[k] = "O"
                continue
            codon = list(model.wt_codon[k])
            for i, off, b in items:
                if cs[i] == "D":
                    codon[off] = b
            codon = "".join(codon)
            if k == tk and model.synonymous_target:
                states[k] = "D" if codon == model.target_design_codon else "W" if codon == model.target_wt_codon else "X"
            else:
                aa = GENETIC_CODE[codon]
                states[k] = "D" if aa == model.design_aa[k] else "W" if aa == model.wt_aa[k] else "X"
        tgt_partial = any(cs[i] == "D" for i, _, _ in tgt) and any(cs[i] == "W" for i, _, _ in tgt)
        if "X" in states.values() or (tgt_partial and states.get(tk) != "D"):
            pX += 1
        elif "O" in states.values():
            pXo += 1
        elif states.get(tk) == "D":
            pD += 1
            if any(cs[i] == "W" for i in model.helpers):
                helpers_wt_in_design = True
        elif states.get(tk) == "W":
            if any(cs[i] == "D" for i in model.helpers):
                pH += 1
            else:
                pW += 1
        elif hasd:
            uH += 1

    if nP:
        call_nt = IE
    elif nD:
        call_nt = MIX if significant_minority(nD, nW, c, alpha) else CE
    elif nO and nO >= nW:
        call_nt = IE
    elif nW >= 2:
        call_nt = NE
    else:
        call_nt = UN
    incomplete = ((call_nt in (CE, IE, MIX) and len(read_d | read_w) < model.n)
                  or (call_nt == NE and len(read_w) < model.n))

    if pX:
        call_target = IE
    elif pXo and pXo >= pW and not pD:
        call_target = IE
    elif pD and pH:
        call_target = MIX
    elif pD:
        if significant_minority(pD, pW, c, alpha):
            call_target = MIX
        else:
            call_target = CE_TARGET if (helpers_wt_in_design or nP) else CE
    elif pH:
        call_target = MIX if significant_minority(pH, pW, c, alpha) else HELPERS_ONLY
    elif pW >= 2:
        call_target = NE
    elif uH and not pW:
        call_target = CE_UNCONFIRMED_TARGET
    else:
        call_target = UN
    return TargetCall(call_nt, call_target, nD, nW, nP, nO, pD, pW, pH, pX, pXo, uH, incomplete,
                      model.synonymous_target)
