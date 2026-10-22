import pytest

from orf_finder.translate import (
    CODON_TABLE,
    STOP_CODONS,
    translate,
    translate_codon,
)


def test_table_covers_every_codon_exactly_once():
    assert len(CODON_TABLE) == 64
    assert len({codon for codon in CODON_TABLE if len(codon) == 3}) == 64


def test_table_has_the_standard_stop_and_start_assignments():
    assert STOP_CODONS == {"TAA", "TAG", "TGA"}
    assert CODON_TABLE["ATG"] == "M"
    assert CODON_TABLE["TGG"] == "W"
    assert sorted(set(CODON_TABLE.values())) == sorted("*ACDEFGHIKLMNPQRSTVWY")


@pytest.mark.parametrize(
    ("codon", "amino_acid"),
    [
        ("TTT", "F"), ("TTA", "L"), ("CTG", "L"), ("ATT", "I"), ("GTG", "V"),
        ("TCA", "S"), ("CCC", "P"), ("ACG", "T"), ("GCC", "A"), ("TAT", "Y"),
        ("CAC", "H"), ("CAG", "Q"), ("AAT", "N"), ("AAA", "K"), ("GAT", "D"),
        ("GAG", "E"), ("TGC", "C"), ("CGA", "R"), ("AGA", "R"), ("AGC", "S"),
        ("GGG", "G"),
    ],
)
def test_representative_codons(codon, amino_acid):
    assert translate_codon(codon) == amino_acid


def test_six_fold_families_are_consistent():
    for codons, amino_acid in (
        (("TTA", "TTG", "CTT", "CTC", "CTA", "CTG"), "L"),
        (("TCT", "TCC", "TCA", "TCG", "AGT", "AGC"), "S"),
        (("CGT", "CGC", "CGA", "CGG", "AGA", "AGG"), "R"),
    ):
        assert {CODON_TABLE[codon] for codon in codons} == {amino_acid}


def test_ambiguous_codons_translate_to_x():
    assert translate_codon("ANG") == "X"
    assert translate_codon("NNN") == "X"


def test_translate_is_case_insensitive_and_drops_a_partial_codon():
    assert translate("atgaaatag") == "MK*"
    assert translate("ATGAAATAGAC") == "MK*"
    assert translate("AT") == ""
