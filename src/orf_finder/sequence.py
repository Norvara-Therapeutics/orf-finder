"""Nucleic-acid normalisation and summary metrics."""

from __future__ import annotations

from dataclasses import dataclass, field

_DNA_COMPLEMENT = str.maketrans(
    "ACGTUNRYSWKMBDHVacgtunryswkmbdhv",
    "TGCAANYRSWMKVHDBtgcaanyrswmkvhdb",
)

#: Residue masses of the nucleoside monophosphates in a polynucleotide chain,
#: as free acids (g/mol). The chain adjustment accounts for a 5'-hydroxyl
#: terminus rather than a 5'-phosphate.
_DNA_RESIDUE_MASS = {"A": 313.21, "C": 289.18, "G": 329.21, "T": 304.20}
_RNA_RESIDUE_MASS = {"A": 329.21, "C": 305.18, "G": 345.21, "U": 306.17}
_TERMINUS_ADJUSTMENT = -61.96
_GAPS = frozenset("-.~")

def reverse_complement(sequence: str) -> str:
    """Reverse-complement a sequence, preserving case and IUPAC ambiguity codes."""
    return sequence.translate(_DNA_COMPLEMENT)[::-1]


def is_rna(sequence: str) -> bool:
    """True when the sequence uses U in place of T."""
    upper = sequence.upper()
    return "U" in upper and "T" not in upper


def normalise(sequence: str) -> str:
    """Upper-case a sequence and transcribe U to T so one codon table serves both."""
    return sequence.upper().replace("U", "T")


@dataclass(frozen=True)
class SequenceMetrics:
    """Composition and mass summary for one sequence.

    ``gc_percent`` is computed over unambiguous bases only (A+C+G+T); it is
    ``None`` when the sequence contains none. ``molecular_weight`` is for the
    single strand as given, ``molecular_weight_duplex`` for the perfect duplex
    and is ``None`` for RNA input. Masses are approximate when ambiguous bases
    are present, each of which is charged the mean residue mass of the four
    standard bases. Alignment gaps count toward length but not mass or ambiguity.
    """

    length: int
    molecule: str
    counts: dict[str, int] = field(default_factory=dict)
    ambiguous: int = 0
    gc_percent: float | None = None
    molecular_weight: float | None = None
    molecular_weight_duplex: float | None = None
    gaps: int = 0

    @property
    def base_order(self) -> tuple[str, ...]:
        """Canonical display order of the counted bases."""
        return ("A", "U", "C", "G") if self.molecule == "RNA" else ("A", "T", "C", "G")


def _single_strand_mass(sequence: str, masses: dict[str, float]) -> float | None:
    mean_residue = sum(masses.values()) / len(masses)
    total = 0.0
    residues = 0
    for base in sequence:
        if base in _GAPS:
            continue
        total += masses.get(base, mean_residue)
        residues += 1
    if not residues:
        return None
    return round(total + _TERMINUS_ADJUSTMENT, 2)


def metrics(sequence: str) -> SequenceMetrics:
    """Compute length, base composition, GC content and estimated mass."""
    upper = sequence.upper()
    molecule = "RNA" if is_rna(upper) else "DNA"
    # A DNA sequence carrying stray U is counted as though those positions were
    # T, which is how find_orfs and translate read them.
    if molecule == "DNA":
        upper = upper.replace("U", "T")
    tracked = ("A", "U", "C", "G") if molecule == "RNA" else ("A", "T", "C", "G")
    counts = {base: upper.count(base) for base in tracked}
    counts["N"] = upper.count("N")
    gaps = sum(upper.count(gap) for gap in _GAPS)
    ambiguous = len(upper) - sum(counts[base] for base in tracked) - gaps

    unambiguous = sum(counts[base] for base in tracked)
    gc = (
        round(100.0 * (counts["C"] + counts["G"]) / unambiguous, 2)
        if unambiguous
        else None
    )

    masses = _RNA_RESIDUE_MASS if molecule == "RNA" else _DNA_RESIDUE_MASS
    weight = _single_strand_mass(upper, masses)
    duplex = None
    if molecule == "DNA" and weight is not None:
        complement = _single_strand_mass(reverse_complement(upper), masses)
        duplex = None if complement is None else round(weight + complement, 2)

    return SequenceMetrics(
        length=len(sequence),
        molecule=molecule,
        counts=counts,
        ambiguous=ambiguous,
        gc_percent=gc,
        molecular_weight=weight,
        molecular_weight_duplex=duplex,
        gaps=gaps,
    )
