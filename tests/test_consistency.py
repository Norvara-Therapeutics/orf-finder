"""Agreement between the modules, and the invariants the scanner relies on."""

import random

import pytest

from orf_finder.orf import _drop_nested, cds_sequence, find_orfs
from orf_finder.sequence import metrics, reverse_complement
from orf_finder.translate import translate, translate_codon

RNA_ORF = "AUGAAACCCGGGUAA"


def test_translate_reads_uracil_like_the_scanner():
    assert translate(RNA_ORF) == "MKPG*"
    assert translate_codon("aug") == "M"
    (orf,) = find_orfs(RNA_ORF, min_length=15, strand="plus")
    assert orf.protein == translate(RNA_ORF)[:-1]


def test_metrics_read_stray_uracil_in_dna_like_the_scanner():
    # Both T and U present, so the molecule is DNA and U is read as T.
    summary = metrics("ATGAAAUAA")
    assert summary.molecule == "DNA"
    assert summary.counts["T"] == 2
    assert summary.ambiguous == 0
    (orf,) = find_orfs("ATGAAAUAA", min_length=9, strand="plus")
    assert orf.protein == "MK"


def test_rna_metrics_count_uracil_natively():
    summary = metrics(RNA_ORF)
    assert summary.molecule == "RNA"
    assert summary.counts["U"] == 2
    assert summary.ambiguous == 0


def test_length_matches_the_reported_span_for_every_orf():
    rng = random.Random(4242)
    sequence = "".join(rng.choice("ACGT") for _ in range(4000))
    for mode in ("atg", "alt", "any"):
        orfs = find_orfs(
            sequence, min_length=30, start_codons=mode, include_truncated=True
        )
        assert orfs
        for orf in orfs:
            assert orf.length_nt == orf.end - orf.start + 1
            assert len(cds_sequence(orf, sequence)) == orf.length_nt


@pytest.mark.parametrize("mode", ["atg", "any"])
def test_the_nesting_sweep_agrees_with_direct_comparison(mode):
    """The sweep replaces an all-pairs scan; they must select the same ORFs."""
    rng = random.Random(99)
    for _ in range(25):
        sequence = "".join(rng.choice("ACGT") for _ in range(rng.randint(90, 900)))
        orfs = find_orfs(sequence, min_length=9, start_codons=mode)
        expected = [
            candidate
            for candidate in orfs
            if not any(
                other.strand == candidate.strand
                and other.length_nt > candidate.length_nt
                and other.start <= candidate.start
                and candidate.end <= other.end
                for other in orfs
            )
        ]
        assert _drop_nested(orfs) == expected


def test_the_nesting_sweep_preserves_order_and_keeps_identical_spans():
    # An ORF on each strand with exactly the same span: neither contains the other.
    sequence = "ATG" + "AAA" * 4 + "TAA"
    both = find_orfs(sequence + reverse_complement(sequence), min_length=18)
    kept = _drop_nested(both)
    assert kept == both
    assert [orf.index for orf in kept] == sorted(orf.index for orf in kept)
