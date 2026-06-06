"""Tests for WGS Stage 3 reconciliation (trace_crispr.wgs.reconcile)."""

import pysam

from trace_crispr.wgs.reconcile import (
    CONTAMINATION,
    LOW_CONF,
    MULTI_EDIT,
    NO_EDIT,
    SINGLE_EDIT,
    SITE_FULL,
    SITE_NOCOV,
    SITE_PARTIAL,
    SITE_WT,
    EditSite,
    SiteCall,
    _designed_offsets,
    build_all_edit_sites,
    build_edit_sites,
    classify_site,
    designed_site_vaf,
    load_oligo_windows,
    reconcile_cell,
    reconciliation_to_row,
    window_edit_vaf,
)


def _sc(status, vaf=None, cov=10, oligo="o", chrom="chrI", start=100, end=110):
    if vaf is None:
        vaf = {SITE_FULL: 1.0, SITE_PARTIAL: 0.5, SITE_WT: 0.0, SITE_NOCOV: 0.0}[status]
    return SiteCall(EditSite(chrom, start, end, oligo), vaf, cov, status)


# --------------------------------------------------------------------------- #
# classify_site thresholds                                                     #
# --------------------------------------------------------------------------- #


def test_classify_site_thresholds():
    assert classify_site(1.0, 10) == SITE_FULL
    assert classify_site(0.50, 10) == SITE_PARTIAL
    assert classify_site(0.00, 10) == SITE_WT
    assert classify_site(1.0, 2) == SITE_NOCOV  # low coverage dominates
    assert classify_site(0.85, 10) == SITE_FULL  # tau_full boundary inclusive
    assert classify_site(0.15, 10) == SITE_WT  # tau_wt boundary inclusive
    # 0.16 at cov 10 == ~2 non-ref reads, below the default read-count floor -> WT
    assert classify_site(0.16, 10) == SITE_WT
    assert classify_site(0.30, 10) == SITE_PARTIAL  # 3 non-ref reads clears the floor


def test_classify_site_read_count_floor():
    # The floor guards the PARTIAL band against single-/few-read false positives.
    assert classify_site(0.25, 4) == SITE_WT  # 1 non-ref read at min_cov -> not PARTIAL
    assert classify_site(0.50, 4) == SITE_WT  # 2 non-ref reads -> still below floor
    assert classify_site(0.30, 60) == SITE_PARTIAL  # genuine secondary (~18 reads)
    # The floor is tunable and does not touch FULL/WT bands.
    assert classify_site(0.25, 4, min_partial_reads=1) == SITE_PARTIAL
    assert classify_site(0.90, 4) == SITE_FULL  # high-VAF call is unaffected by the floor


# --------------------------------------------------------------------------- #
# _designed_offsets — substitution / deletion / insertion edge cases          #
# --------------------------------------------------------------------------- #


def test_designed_offsets_substitution():
    assert _designed_offsets("ACGT", "ATGT") == [1]  # single sub
    assert _designed_offsets("ACGT", "ATTT") == [1, 2, 3]  # multi-base sub block
    assert _designed_offsets("ACGT", "ACGT") == []  # no change


def test_designed_offsets_deletion():
    assert _designed_offsets("ACGT", "AGT") == [1]  # delete C (offset 1)
    assert _designed_offsets("ACGT", "AT") == [1, 2]  # delete CG


def test_designed_offsets_insertion_anchors_to_prior_base():
    # Insertion evidence is read from PileupRead.indel at the anchor (prior) base.
    assert _designed_offsets("ACGT", "ACXGT") == [1]  # insert between C and G -> anchor 1
    assert _designed_offsets("ACGT", "XACGT") == [0]  # insert at start -> clamp to 0


# --------------------------------------------------------------------------- #
# reconcile_cell — every label branch                                          #
# --------------------------------------------------------------------------- #


def test_single_edit():
    c = reconcile_cell("s", [_sc(SITE_FULL, start=100)])
    assert c.label == SINGLE_EDIT and c.confidence == "high" and c.n_full == 1
    assert c.n_sites == 1 and c.n_covered == 1


def test_multi_edit():
    c = reconcile_cell("s", [_sc(SITE_FULL, start=100), _sc(SITE_FULL, start=200, oligo="o2")])
    assert c.label == MULTI_EDIT and c.n_full == 2 and c.confidence == "high"


def test_partial_vaf_is_contamination_and_overrides_multi():
    sites = [
        _sc(SITE_FULL, start=100),
        _sc(SITE_FULL, start=200, oligo="o2"),
        _sc(SITE_PARTIAL, start=300, oligo="o3"),
    ]
    c = reconcile_cell("s", sites)
    assert c.label == CONTAMINATION and "PARTIAL_VAF" in c.flags


def test_genome_multiplicity_is_contamination():
    # two clean full edits but >=2 genome-track bc1 -> mixture (only one integrates)
    c = reconcile_cell("s", [_sc(SITE_FULL, start=100)], n_genome_bc=2)
    assert c.label == CONTAMINATION and "MULTI_GENOME_BC" in c.flags


def test_no_edit():
    c = reconcile_cell("s", [_sc(SITE_WT, start=100)])
    assert c.label == NO_EDIT and c.confidence == "high"
    assert c.n_wt == 1


def test_no_edit_with_uncovered_site_is_low_confidence():
    c = reconcile_cell("s", [_sc(SITE_WT, start=100), _sc(SITE_NOCOV, start=200, oligo="o2")])
    assert c.label == NO_EDIT and c.confidence == "low" and "UNCOVERED_SITES" in c.flags


def test_nothing_covered_is_low_conf():
    c = reconcile_cell("s", [_sc(SITE_NOCOV, start=100)])
    assert c.label == LOW_CONF and c.confidence == "low"


def test_low_bc_confidence_downgrades_a_clean_call():
    c = reconcile_cell("s", [_sc(SITE_FULL, start=100)], bc_confidence="low")
    assert c.label == SINGLE_EDIT and c.confidence == "low"


def test_row_roundtrips_key_fields():
    c = reconcile_cell(
        "s",
        [
            _sc(SITE_FULL, start=100, oligo="oA"),
            _sc(SITE_WT, start=200, end=210, oligo="oB"),
        ],
        n_plasmid_bc=1,
    )
    row = reconciliation_to_row(c)
    assert row[0] == "s" and row[1] == SINGLE_EDIT
    assert row[5:11] == ["2", "2", "1", "0", "1", "0"]
    assert "oA@chrI:100-110=1.00/10:full" in row[11]
    assert "oB@chrI:200-210=0.00/10:wt" not in row[11]
    assert "oB@chrI:200-210=0.00/10:wt" in reconciliation_to_row(c, emit_all_sites=True)[11]


# --------------------------------------------------------------------------- #
# build_edit_sites + load_oligo_windows                                        #
# --------------------------------------------------------------------------- #


def test_build_edit_sites_dedups_by_window():
    ow = {
        "oA": EditSite("chrI", 100, 110, "oA"),
        "oA_dup": EditSite("chrI", 100, 110, "oA_dup"),
        "oB": EditSite("chrII", 200, 210, "oB"),
    }
    sites = build_edit_sites(["oA", "oA", "oB", "missing"], ow)
    assert len(sites) == 2
    assert {s.chrom for s in sites} == {"chrI", "chrII"}


def test_build_all_edit_sites_preserves_distinct_designs_in_same_window():
    ow = {
        "oA": EditSite("chrI", 100, 110, "oA", wt_seq="AAAAAAAAAAA", design_donor="AAAAGAAAAAA"),
        "oB": EditSite("chrI", 100, 110, "oB", wt_seq="AAAAAAAAAAA", design_donor="AAAAATAAAAA"),
        "oA_dup": EditSite(
            "chrI", 100, 110, "oA_dup", wt_seq="AAAAAAAAAAA", design_donor="AAAAGAAAAAA"
        ),
    }
    sites = build_all_edit_sites(ow)
    assert len(sites) == 2


def test_load_oligo_windows(tmp_path):
    p = tmp_path / "win.tsv"
    p.write_text(
        "oligo_name\tchrom\tdonor_start_coord\tdonor_end_coord\twt_seq\tdesign_donor\n"
        "oA\tchrI\t100\t110\tAAAAAAAAAAA\tAAAAGAAAAAA\n"
        "oB\tchrII\t200\t210\tCCCCCCCCCCC\tCCCCGCCCCCC\n"
        "bad\tchrIII\tNA\t300\tAAAA\tAAAT\n"  # non-numeric -> skipped
    )
    ow = load_oligo_windows(str(p))  # synthetic_donor absent -> falls back to design_donor
    assert ow["oA"] == EditSite(
        "chrI", 100, 110, "oA", wt_seq="AAAAAAAAAAA", design_donor="AAAAGAAAAAA"
    )
    assert ow["oB"].chrom == "chrII"


def test_load_oligo_windows_prefers_synthetic_donor(tmp_path):
    # synthetic_donor is the integrated edit; design_donor is silent (== wt) for ~90%
    # of oligos, so when both are present synthetic_donor must win.
    p = tmp_path / "win.tsv"
    p.write_text(
        "oligo_name\tchrom\tdonor_start_coord\tdonor_end_coord\twt_seq\tdesign_donor\tsynthetic_donor\n"
        # design_donor identical to wt_seq (silent); synthetic_donor carries the real edit
        "oA\tchrI\t100\t110\tAAAAAAAAAAA\tAAAAAAAAAAA\tAAAAGTAAAAA\n"
    )
    ow = load_oligo_windows(str(p))
    assert ow["oA"].design_donor == "AAAAGTAAAAA"
    # explicit override still works
    ow2 = load_oligo_windows(str(p), design_col="design_donor")
    assert ow2["oA"].design_donor == "AAAAAAAAAAA"


# --------------------------------------------------------------------------- #
# window_edit_vaf on a crafted BAM                                             #
# --------------------------------------------------------------------------- #


def _bam(tmp_path, refseq, reads):
    hdr = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chrI", "LN": len(refseq)}]}
    )
    u = tmp_path / "u.bam"
    with pysam.AlignmentFile(str(u), "wb", header=hdr) as out:
        for i, (pos1, seq) in enumerate(reads):
            a = pysam.AlignedSegment(hdr)
            a.query_name = f"r{i}"
            a.flag = 0
            a.reference_id = 0
            a.reference_start = pos1 - 1
            a.mapping_quality = 60
            a.cigarstring = f"{len(seq)}M"
            a.query_sequence = seq
            a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
            out.write(a)
    s = tmp_path / "s.bam"
    pysam.sort("-o", str(s), str(u))
    pysam.index(str(s))
    return str(s)


def test_window_edit_vaf_half_edited(tmp_path):
    ref = "A" * 40
    fa = tmp_path / "ref.fa"
    fa.write_text(">chrI\n" + ref + "\n")
    pysam.faidx(str(fa))
    wt = "A" * 20
    mut = "A" * 9 + "G" + "A" * 10  # mismatch at 1-based pos 10
    bam_path = _bam(tmp_path, ref, [(1, wt)] * 4 + [(1, mut)] * 4)
    with pysam.AlignmentFile(bam_path, "rb") as bam, pysam.FastaFile(str(fa)) as fasta:
        vaf, cov = window_edit_vaf(bam, fasta, "chrI", 5, 15)
    assert abs(vaf - 0.5) < 1e-9
    assert cov == 8


def test_window_edit_vaf_wt_is_zero(tmp_path):
    ref = "A" * 40
    fa = tmp_path / "ref.fa"
    fa.write_text(">chrI\n" + ref + "\n")
    pysam.faidx(str(fa))
    bam_path = _bam(tmp_path, ref, [(1, "A" * 20)] * 6)
    with pysam.AlignmentFile(bam_path, "rb") as bam, pysam.FastaFile(str(fa)) as fasta:
        vaf, cov = window_edit_vaf(bam, fasta, "chrI", 5, 15)
    assert vaf == 0.0 and cov == 6


def test_designed_site_vaf_ignores_background_mismatch(tmp_path):
    ref = "A" * 40
    fa = tmp_path / "ref.fa"
    fa.write_text(">chrI\n" + ref + "\n")
    pysam.faidx(str(fa))
    background_only = "A" * 5 + "C" + "A" * 14  # non-designed mismatch at 1-based pos 6
    bam_path = _bam(tmp_path, ref, [(1, background_only)] * 6)
    site = EditSite(
        "chrI",
        5,
        15,
        "oA",
        wt_seq="A" * 11,
        design_donor="A" * 5 + "G" + "A" * 5,  # designed mismatch at 1-based pos 10
    )
    with pysam.AlignmentFile(bam_path, "rb") as bam, pysam.FastaFile(str(fa)) as fasta:
        designed_vaf, designed_cov = designed_site_vaf(bam, fasta, site)
        window_vaf, window_cov = window_edit_vaf(bam, fasta, "chrI", 5, 15)
    assert designed_vaf == 0.0 and designed_cov == 6
    assert window_vaf == 1.0 and window_cov == 6


def test_designed_site_vaf_orients_minus_strand_design(tmp_path):
    ref = "T" * 4 + "AGCC" + "T" * 4
    edited = "T" * 4 + "AGCT" + "T" * 4
    fa = tmp_path / "ref.fa"
    fa.write_text(">chrI\n" + ref + "\n")
    pysam.faidx(str(fa))
    bam_path = _bam(tmp_path, ref, [(1, edited)] * 6)
    site = EditSite(
        "chrI",
        5,
        8,
        "oMinus",
        wt_seq="GGCT",       # reverse-complement matches AGCC reference window
        design_donor="AGCT",  # reverse-complement changes ref-window last base C -> T
    )
    with pysam.AlignmentFile(bam_path, "rb") as bam, pysam.FastaFile(str(fa)) as fasta:
        vaf, cov = designed_site_vaf(bam, fasta, site)
    assert vaf == 1.0 and cov == 6


def test_designed_site_vaf_tolerates_strain_snp_in_orientation(tmp_path):
    # Real libraries: wt_seq is reverse-complement of the reference with a few strain
    # SNPs / synthetic recoding (not an exact match). Orientation must still resolve and
    # the designed offset must be scored independent of the background SNP.
    ref = "C" * 4 + "TTTAAGCC" + "C" * 4  # window chrI:5-12 == "TTTAAGCC" (strain SNP @ off3)
    edited = "TTTATGCC"  # designed edit at window offset 4 (A->T); off3 stays ref
    fa = tmp_path / "ref.fa"
    fa.write_text(">chrI\n" + ref + "\n")
    pysam.faidx(str(fa))
    bam_path = _bam(tmp_path, ref, [(5, edited)] * 6)
    site = EditSite(
        "chrI", 5, 12, "oRCsnp",
        wt_seq="GGCTAAAA",      # rc == "TTTTAGCC": 1 mismatch vs ref (the strain SNP)
        design_donor="GGCAAAAA",  # rc edit changes window offset 4 only
    )
    with pysam.AlignmentFile(bam_path, "rb") as bam, pysam.FastaFile(str(fa)) as fasta:
        vaf, cov = designed_site_vaf(bam, fasta, site)
    assert vaf == 1.0 and cov == 6


def test_site_design_offsets_raises_on_wrong_coords(tmp_path):
    from trace_crispr.wgs.reconcile import _site_design_offsets

    ref = "ACGT" * 10
    fa = tmp_path / "ref.fa"
    fa.write_text(">chrI\n" + ref + "\n")
    pysam.faidx(str(fa))
    site = EditSite(
        "chrI", 5, 12, "oBad",
        wt_seq="TTTTTTTT",  # matches neither orientation of the ACGT-repeat window
        design_donor="TTTTTTTA",
    )
    with pysam.FastaFile(str(fa)) as fasta:
        try:
            _site_design_offsets(fasta, site)
        except ValueError:
            pass
        else:
            raise AssertionError("expected ValueError on unmappable wt_seq")
