"""Codon translation against any NCBI genetic code."""

from __future__ import annotations

from orf_finder.genetic_code import (
    CODONS,
    GeneticCode,
    get_genetic_code,
)
from orf_finder.sequence import normalise

#: The standard genetic code, used whenever no table is named.
STANDARD_CODE: GeneticCode = get_genetic_code(1)

#: Codon to amino acid under the standard code, kept for convenience.
CODON_TABLE = dict(STANDARD_CODE.codons)

#: Stop codons of the standard code.
STOP_CODONS = STANDARD_CODE.stop_codons

UNKNOWN_AMINO_ACID = "X"

__all__ = [
    "CODONS",
    "CODON_TABLE",
    "STANDARD_CODE",
    "STOP_CODONS",
    "UNKNOWN_AMINO_ACID",
    "translate",
    "translate_codon",
]


def translate_codon(codon: str, code: GeneticCode = STANDARD_CODE) -> str:
    """Translate one codon, returning 'X' for anything containing an ambiguity code.

    Case is ignored and U is read as T, so RNA codons translate as their DNA
    equivalents.
    """
    return code.codons.get(normalise(codon), UNKNOWN_AMINO_ACID)


def translate(sequence: str, code: GeneticCode = STANDARD_CODE) -> str:
    """Translate a nucleotide sequence in frame 1, ignoring a trailing partial codon.

    Case is ignored and U is read as T, matching :func:`orf_finder.orf.find_orfs`.
    """
    seq = normalise(sequence)
    codons = code.codons
    return "".join(
        codons.get(seq[i : i + 3], UNKNOWN_AMINO_ACID)
        for i in range(0, len(seq) - len(seq) % 3, 3)
    )
