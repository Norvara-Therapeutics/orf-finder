"""Start-codon modes, nesting, fragments and coding-sequence extraction."""

import pytest

from orf_finder.orf import cds_sequence, find_orfs
from orf_finder.sequence import reverse_complement

# Frame +1 holds two stop-to-stop regions; frame +1 start codon at base 10.
SEQUENCE = "AAACCCTAA" + "ATGAAACCCGGG" + "TAA" + "CCCAAAGGG"


def test_atg_mode_requires_atg():
    (orf,) = find_orfs(SEQUENCE, min_length=15, start_codons="atg", strand="plus")
    assert (orf.start, orf.end) == (10, 24)
    assert orf.protein == "MKPG"
    assert orf.five_prime_partial is False


def test_alt_mode_adds_the_tables_initiation_codons():
    sequence = "CTG" + "AAACCCGGG" + "TAA"
    assert find_orfs(sequence, min_length=15, start_codons="atg", strand="plus") == []
    (orf,) = find_orfs(sequence, min_length=15, start_codons="alt", strand="plus")
    assert orf.start_codon == "CTG"
    assert orf.protein == "MKPG", "an initiation codon translates as methionine"


def test_an_explicit_codon_list_still_works():
    sequence = "GTG" + "AAACCCGGG" + "TAA"
    (orf,) = find_orfs(
        sequence, min_length=15, start_codons=["gtg"], strand="plus"
    )
    assert orf.start_codon == "GTG"
    # GTG does not initiate under the standard code, so it keeps its own residue.
    assert orf.protein == "VKPG"
    assert orf.has_initiation_codon is False
    assert orf.complete is False


def test_any_sense_mode_reports_the_regions_between_stops():
    orfs = find_orfs(SEQUENCE, min_length=3, start_codons="any", strand="plus")
    frame_one = sorted((o for o in orfs if o.frame == 1), key=lambda o: o.start)
    assert [(o.start, o.end) for o in frame_one] == [(1, 9), (10, 24)]
    # The first region of a frame has no upstream stop.
    assert frame_one[0].five_prime_partial is True
    assert frame_one[0].complete is False
    assert frame_one[0].protein == "KP"
    assert frame_one[1].five_prime_partial is False
    # In any-sense mode the first codon is not an initiator.
    assert frame_one[1].protein == "MKPG"


def test_any_sense_mode_does_not_need_a_start_codon():
    sequence = "AAACCCTAA"
    assert find_orfs(sequence, min_length=3, start_codons="atg", strand="plus") == []
    (orf,) = find_orfs(sequence, min_length=3, start_codons="any", strand="plus")
    assert orf.protein == "KP"
    assert orf.stop_codon == "TAA"


def test_a_stop_bounded_any_sense_region_is_not_initiated():
    (orf,) = find_orfs("TAACCCTAA", min_length=6, start_codons="any", strand="plus")
    assert orf.start_codon == "CCC"
    assert orf.five_prime_partial is False
    assert orf.has_initiation_codon is False
    assert orf.complete is False


def test_ignore_nested_drops_orfs_inside_a_longer_same_strand_orf():
    # A 300 nt ORF in frame +1 with a shorter ORF wholly inside it in frame +2.
    outer = "ATG" + "AAA" * 98 + "TAA"
    inner_start = 4  # frame +2
    sequence = outer[:inner_start] + "ATG" + "CCC" * 20 + "TGA" + outer[inner_start + 66:]
    orfs = find_orfs(sequence, min_length=60, strand="plus")
    assert len(orfs) >= 2
    kept = find_orfs(sequence, min_length=60, strand="plus", ignore_nested=True)
    assert len(kept) == 1
    assert kept[0].length_nt == max(orf.length_nt for orf in orfs)
    assert kept[0].label == "ORF1"


def test_ignore_nested_is_judged_within_a_strand():
    plus = "ATG" + "AAA" * 20 + "TAA"
    minus = reverse_complement("ATG" + "CCC" * 5 + "TAA")
    sequence = plus[:6] + minus + plus[6 + len(minus):]
    orfs = find_orfs(sequence, min_length=18, ignore_nested=True)
    strands = {orf.strand for orf in orfs}
    assert strands == {"+", "-"}, "a minus-strand ORF inside a plus-strand ORF is kept"


def test_offset_shifts_reported_coordinates():
    (plain,) = find_orfs(SEQUENCE, min_length=15, strand="plus")
    (shifted,) = find_orfs(SEQUENCE, min_length=15, strand="plus", offset=1000)
    assert (shifted.start, shifted.end) == (plain.start + 1000, plain.end + 1000)


def test_cds_sequence_returns_the_coding_strand():
    (orf,) = find_orfs(SEQUENCE, min_length=15, strand="plus")
    assert cds_sequence(orf, SEQUENCE) == "ATGAAACCCGGGTAA"
    reverse = reverse_complement(SEQUENCE)
    (minus,) = find_orfs(reverse, min_length=15, strand="minus")
    assert cds_sequence(minus, reverse) == "ATGAAACCCGGGTAA"


def test_cds_sequence_honours_the_offset():
    (orf,) = find_orfs(SEQUENCE, min_length=15, strand="plus", offset=50)
    assert cds_sequence(orf, SEQUENCE, offset=50) == "ATGAAACCCGGGTAA"


@pytest.mark.parametrize("mode", ["ATG", "Any", "alt"])
def test_start_modes_are_case_insensitive(mode):
    find_orfs(SEQUENCE, min_length=3, start_codons=mode)


def test_an_unknown_start_mode_is_rejected():
    with pytest.raises(ValueError, match="start mode must be one of"):
        find_orfs(SEQUENCE, start_codons="sometimes")


def test_an_unknown_table_is_rejected():
    with pytest.raises(KeyError):
        find_orfs(SEQUENCE, table=7)
