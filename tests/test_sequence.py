import pytest

from orf_finder.sequence import (
    is_rna,
    metrics,
    normalise,
    reverse_complement,
)


def test_reverse_complement_handles_case_and_ambiguity():
    assert reverse_complement("ACGT") == "ACGT"
    assert reverse_complement("AAAC") == "GTTT"
    assert reverse_complement("acgt") == "acgt"
    assert reverse_complement("ACGTN") == "NACGT"
    assert reverse_complement("ACGTRYKM") == "KMRYACGT"
    assert reverse_complement("AUGC") == "GCAT"


def test_rna_detection_and_normalisation():
    assert is_rna("AUGC")
    assert not is_rna("ATGC")
    assert not is_rna("")
    assert not is_rna("ATGCU")  # both T and U present: treated as DNA
    assert normalise("augc") == "ATGC"


def test_counts_gc_and_length():
    summary = metrics("AACCGGTTNN")
    assert summary.length == 10
    assert summary.molecule == "DNA"
    assert summary.counts == {"A": 2, "T": 2, "C": 2, "G": 2, "N": 2}
    assert summary.ambiguous == 2
    assert summary.gc_percent == pytest.approx(50.0)
    assert summary.base_order == ("A", "T", "C", "G")


def test_gc_percent_ignores_ambiguous_bases():
    assert metrics("GGNN").gc_percent == pytest.approx(100.0)
    assert metrics("NNNN").gc_percent is None


def test_molecular_weight_follows_the_documented_residue_masses():
    # 313.21 + 304.20 + 289.18 + 329.21 - 61.96
    assert metrics("ATGC").molecular_weight == pytest.approx(1173.84)
    # The duplex of a self-complementary sequence is twice the single strand.
    assert metrics("ATGC").molecular_weight_duplex == pytest.approx(2347.68)


def test_rna_metrics_use_uracil_and_report_no_duplex_mass():
    summary = metrics("AUGC")
    assert summary.molecule == "RNA"
    assert summary.base_order == ("A", "U", "C", "G")
    assert summary.counts["U"] == 1
    assert "T" not in summary.counts
    # 329.21 + 306.17 + 345.21 + 305.18 - 61.96
    assert summary.molecular_weight == pytest.approx(1223.81)
    assert summary.molecular_weight_duplex is None


def test_empty_sequence_reports_no_mass():
    summary = metrics("")
    assert summary.length == 0
    assert summary.molecular_weight is None
    assert summary.gc_percent is None


def test_alignment_gaps_do_not_add_residue_mass_or_ambiguity():
    plain = metrics("ACGT")
    gapped = metrics("A-.~CGT")
    assert gapped.length == 7
    assert gapped.gaps == 3
    assert gapped.ambiguous == 0
    assert gapped.molecular_weight == plain.molecular_weight
    assert gapped.molecular_weight_duplex == plain.molecular_weight_duplex
    assert metrics("-.~").molecular_weight is None
