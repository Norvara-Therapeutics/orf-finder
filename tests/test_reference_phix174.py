"""Checks against an annotated reference.

Gene J of Escherichia phage phiX174 (RefSeq NC_001422.1, nucleotides 848..964)
and its recorded protein product. The translation of every annotated coding
sequence in that record reproduces the deposited protein exactly; gene J is
committed here as a compact offline check of the genetic code, the frame
arithmetic and the strand mapping.
"""

from orf_finder.orf import find_orfs
from orf_finder.sequence import reverse_complement
from orf_finder.translate import translate

GENE_J = "ATGTCTAAAGGTAAAAAACGTTCTGGCGCTCGCCCTGGTCGTCCGCAGCCGTTGCGAGGTACTAAAGGCAAGCGTAAAGGCGCTCGTCTTTGGTATGTAGGTGGTCAACAATTTTAA"
GENE_J_PROTEIN = "MSKGKKRSGARPGRPQPLRGTKGKRKGARLWYVGGQQF"


def test_translation_reproduces_the_deposited_protein():
    assert translate(GENE_J) == GENE_J_PROTEIN + "*"


def test_the_gene_is_recovered_as_a_plus_strand_orf():
    (orf,) = find_orfs(GENE_J, min_length=len(GENE_J), strand="plus")
    assert (orf.start, orf.end, orf.length_nt) == (1, 117, 117)
    assert (orf.frame, orf.strand) == (1, "+")
    assert orf.protein == GENE_J_PROTEIN
    assert orf.length_aa == 38
    assert (orf.start_codon, orf.stop_codon) == ("ATG", "TAA")


def test_the_gene_is_recovered_on_the_minus_strand_with_mirrored_coordinates():
    (orf,) = find_orfs(reverse_complement(GENE_J), min_length=len(GENE_J), strand="minus")
    assert (orf.start, orf.end) == (1, 117)
    assert orf.frame == -1
    assert orf.protein == GENE_J_PROTEIN
