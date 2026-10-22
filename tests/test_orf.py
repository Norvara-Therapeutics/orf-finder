import random

import pytest

from orf_finder.orf import find_orfs
from orf_finder.sequence import reverse_complement

PLUS_ORF = "ATG" + "AAA" * 4 + "TAA"  # 18 nt, peptide MKKKK
FLANKED = "CC" + PLUS_ORF + "GG"


@pytest.mark.parametrize("table,stop", [(27, "TGA"), (28, "TAA"), (31, "TAG")])
def test_context_dependent_terminal_codon_is_excluded_from_peptide(table, stop):
    (orf,) = find_orfs("ATGAAA" + stop, min_length=9, strand="plus", table=table)
    assert orf.stop_codon == stop
    assert orf.protein == "MK"
    assert orf.length_aa == 2


def test_finds_a_plus_strand_orf_with_inclusive_coordinates():
    (orf,) = find_orfs(FLANKED, min_length=18)
    assert (orf.strand, orf.frame) == ("+", 3)
    assert (orf.start, orf.end) == (3, 20)
    assert orf.length_nt == 18
    assert orf.length_aa == 5
    assert orf.protein == "MKKKK"
    assert (orf.start_codon, orf.stop_codon) == ("ATG", "TAA")
    assert orf.truncated is False
    assert orf.label == "ORF1"


def test_minus_strand_coordinates_mirror_the_plus_strand():
    (orf,) = find_orfs(reverse_complement(FLANKED), min_length=18)
    assert (orf.strand, orf.frame) == ("-", -3)
    assert (orf.start, orf.end) == (3, 20)
    assert orf.protein == "MKKKK"


def test_min_length_counts_the_stop_codon():
    assert len(find_orfs(FLANKED, min_length=18)) == 1
    assert find_orfs(FLANKED, min_length=19) == []


def test_one_orf_is_reported_per_stop_from_the_first_start():
    (orf,) = find_orfs("ATGATGAAATAA", min_length=3)
    assert orf.start == 1
    assert orf.protein == "MMK"


def test_a_new_orf_may_open_after_a_stop_in_the_same_frame():
    orfs = find_orfs(PLUS_ORF + PLUS_ORF, min_length=18)
    assert [(orf.start, orf.end) for orf in orfs] == [(1, 18), (19, 36)]


def test_truncated_orfs_are_excluded_unless_requested():
    assert find_orfs("ATGAAAAAA", min_length=3) == []
    (orf,) = find_orfs("ATGAAAAAA", min_length=3, include_truncated=True)
    assert orf.truncated is True
    assert orf.stop_codon is None
    assert orf.protein == "MKK"
    assert (orf.start, orf.end, orf.length_nt) == (1, 9, 9)


def test_a_truncated_orf_stops_at_the_last_whole_codon():
    (orf,) = find_orfs("ATGAAAAAAG", min_length=3, include_truncated=True)
    assert (orf.end, orf.length_nt) == (9, 9)


def test_strand_selection():
    sequence = reverse_complement(FLANKED)
    assert find_orfs(sequence, min_length=18, strand="plus") == []
    assert len(find_orfs(sequence, min_length=18, strand="minus")) == 1


def test_alternative_start_codons_are_honoured():
    sequence = "GTG" + "AAA" * 4 + "TAA"
    assert find_orfs(sequence, min_length=18) == []
    (orf,) = find_orfs(sequence, min_length=18, start_codons=("ATG", "GTG"))
    assert orf.start_codon == "GTG"
    assert orf.protein == "VKKKK"


def test_uracil_and_lower_case_input_are_accepted():
    (orf,) = find_orfs(PLUS_ORF.replace("T", "U").lower(), min_length=18)
    assert orf.protein == "MKKKK"


def test_orfs_are_sorted_by_descending_length_and_numbered_in_order():
    short = "ATG" + "AAA" + "TAA"
    sequence = short + "C" + PLUS_ORF
    orfs = find_orfs(sequence, min_length=6)
    assert [orf.length_nt for orf in orfs] == [18, 9]
    assert [orf.index for orf in orfs] == [1, 2]
    assert [orf.label for orf in orfs] == ["ORF1", "ORF2"]


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [
        ({"strand": "either"}, "strand must be one of"),
        ({"min_length": -1}, "must not be negative"),
        ({"start_codons": ()}, "at least one start codon"),
        ({"start_codons": ("AT",)}, "not a codon"),
        ({"start_codons": ("ANG",)}, "not a codon"),
    ],
)
def test_invalid_arguments_are_rejected(kwargs, match):
    with pytest.raises(ValueError, match=match):
        find_orfs(FLANKED, **kwargs)


def test_reported_coordinates_round_trip_to_the_reported_peptide():
    rng = random.Random(20261021)
    sequence = "".join(rng.choice("ACGT") for _ in range(6000))
    orfs = find_orfs(sequence, min_length=60, include_truncated=True)
    assert orfs, "the fixture sequence should contain ORFs"
    for orf in orfs:
        span = sequence[orf.start - 1 : orf.end]
        if orf.strand == "-":
            span = reverse_complement(span)
        assert len(span) == orf.length_nt == orf.end - orf.start + 1
        assert span[:3] == orf.start_codon
        if orf.truncated:
            assert len(orf.protein) * 3 == orf.length_nt
        else:
            assert span[-3:] == orf.stop_codon
            assert len(orf.protein) * 3 == orf.length_nt - 3
        if orf.strand == "+":
            assert (orf.start - 1) % 3 == orf.frame - 1
