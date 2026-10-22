"""Six-frame open reading frame detection.

Frames are numbered +1, +2, +3 by offset into the given strand and -1, -2, -3
by offset into its reverse complement.

Three start-codon modes are supported:

``atg``
    ATG only.
``alt``
    ATG together with every alternative initiation codon the chosen translation
    table declares.
``any``
    Any sense codon, which makes every region between two stop codons an ORF.
    The first region of each frame has no upstream stop and is therefore
    5'-incomplete. In this mode ``start_codon`` holds the first codon of the
    region, which carries no initiation meaning.

In the two start-codon modes one ORF is reported per in-frame stop codon,
beginning at the 5'-most start codon that follows the preceding stop, so nested
downstream starts sharing a stop are not reported separately. An initiation
codon is translated as methionine, as every NCBI table marks its starts and as
deposited coding-sequence translations read; in ``any`` mode the first codon is
not an initiator and is translated normally.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Iterable, Iterator, Sequence

from orf_finder.genetic_code import GeneticCode, get_genetic_code
from orf_finder.sequence import normalise, reverse_complement
from orf_finder.translate import STANDARD_CODE, translate

STRANDS = ("both", "plus", "minus")

#: Start-codon modes accepted by :func:`find_orfs`.
START_MODES = ("atg", "alt", "any")


@dataclass(frozen=True)
class ORF:
    """One open reading frame.

    Coordinates are 1-based and inclusive on the forward strand of the input,
    with ``start`` <= ``end`` on both strands; a minus-strand ORF is read from
    ``end`` towards ``start``. ``length_nt`` spans the whole ORF including its
    stop codon, and ``protein`` excludes the stop. ``has_initiation_codon``
    records whether the selected first codon is an initiator in the chosen
    genetic code; it is false for every ORF found in ``any`` mode.
    """

    index: int
    strand: str
    frame: int
    start: int
    end: int
    length_nt: int
    start_codon: str
    stop_codon: str | None
    protein: str
    truncated: bool
    five_prime_partial: bool = False
    has_initiation_codon: bool = True

    @property
    def length_aa(self) -> int:
        """Length of the translated peptide, excluding any stop codon."""
        return len(self.protein)

    @property
    def label(self) -> str:
        """Short identifier of the form 'ORF1'."""
        return f"ORF{self.index}"

    @property
    def complete(self) -> bool:
        """True when the ORF has both an initiation codon and a stop codon."""
        return (
            self.has_initiation_codon
            and not self.truncated
            and not self.five_prime_partial
        )

    def as_dict(self) -> dict[str, object]:
        """A flat, JSON-serialisable view."""
        return {
            "label": self.label,
            "strand": self.strand,
            "frame": self.frame,
            "start": self.start,
            "end": self.end,
            "length_nt": self.length_nt,
            "length_aa": self.length_aa,
            "start_codon": self.start_codon,
            "stop_codon": self.stop_codon,
            "truncated": self.truncated,
            "five_prime_partial": self.five_prime_partial,
            "has_initiation_codon": self.has_initiation_codon,
            "protein": self.protein,
        }


def cds_sequence(orf: ORF, sequence: str, offset: int = 0) -> str:
    """Return the coding sequence of ``orf``, read 5' to 3' on its own strand.

    Args:
        orf: an ORF reported for ``sequence``.
        sequence: the sequence the ORF was found in; returned bases are
            normalised to upper-case DNA so that they translate to the
            reported peptide.
        offset: the value of ``offset`` passed to :func:`find_orfs`.
    """
    span = normalise(sequence)[orf.start - 1 - offset : orf.end - offset]
    return reverse_complement(span) if orf.strand == "-" else span


def _hit(
    begin: int,
    finish: int,
    total: int,
    sign: int,
    work: str,
    stop_codon: str | None,
    five_prime_partial: bool,
) -> tuple[int, int, str, str | None, str, bool, bool]:
    start, end = (
        (begin + 1, finish) if sign > 0 else (total - finish + 1, total - begin)
    )
    return (
        start,
        end,
        work[begin : begin + 3],
        stop_codon,
        work[begin:finish],
        stop_codon is None,
        five_prime_partial,
    )


def _scan_strand(
    work: str,
    total: int,
    sign: int,
    starts: frozenset[str],
    stops: frozenset[str],
    any_sense: bool,
    include_truncated: bool,
) -> Iterator[tuple[int, tuple[int, int, str, str | None, str, bool, bool]]]:
    """Yield (frame, hit) pairs for one strand."""
    length = len(work)
    for offset in range(3):
        frame = sign * (offset + 1)
        open_at: int | None = offset if any_sense else None
        partial = any_sense  # the first region of the frame has no upstream stop
        for i in range(offset, length - 2, 3):
            codon = work[i : i + 3]
            if codon in stops:
                if open_at is not None and i > open_at:
                    yield frame, _hit(open_at, i + 3, total, sign, work, codon, partial)
                if any_sense:
                    open_at, partial = i + 3, False
                else:
                    open_at = None
            elif open_at is None and codon in starts:
                open_at = i
        if open_at is not None and include_truncated:
            stop_at = offset + 3 * ((length - offset) // 3)
            if stop_at > open_at:
                yield frame, _hit(open_at, stop_at, total, sign, work, None, partial)


def _drop_nested(orfs: list[ORF]) -> list[ORF]:
    """Remove every ORF completely placed within a longer ORF on the same strand.

    Because ``length_nt`` equals ``end - start + 1``, an ORF that spans another
    is automatically at least as long, and is strictly longer unless the two
    spans are identical. Containment is therefore decided by a sweep over each
    strand in start order: an ORF is dropped when some earlier-starting ORF
    reaches at least as far, or when another ORF with the same start reaches
    further. The input order is preserved in the result.
    """
    contained: set[int] = set()
    for strand in ("+", "-"):
        indexed = [(orf, i) for i, orf in enumerate(orfs) if orf.strand == strand]
        indexed.sort(key=lambda pair: (pair[0].start, -pair[0].end))
        furthest_earlier = 0  # furthest end among ORFs with a smaller start
        group_start: int | None = None
        group_end = 0  # furthest end among ORFs sharing the current start
        for orf, position in indexed:
            if orf.start != group_start:
                furthest_earlier = max(furthest_earlier, group_end)
                group_start, group_end = orf.start, orf.end
            if furthest_earlier >= orf.end or group_end > orf.end:
                contained.add(position)
    return [orf for i, orf in enumerate(orfs) if i not in contained]


def _resolve_starts(
    start_codons: Sequence[str] | Iterable[str] | str | None,
    code: GeneticCode,
) -> tuple[frozenset[str], bool]:
    """Resolve the start-codon specification to a codon set and the any-sense flag."""
    spec = "atg" if start_codons is None else start_codons
    if isinstance(spec, str):
        mode = spec.strip().lower()
        if mode not in START_MODES:
            raise ValueError(
                f"start mode must be one of {', '.join(START_MODES)}, "
                f"or a list of codons; got {spec!r}"
            )
        if mode == "any":
            return frozenset(), True
        if mode == "atg":
            return frozenset({"ATG"}), False
        return frozenset({"ATG"}) | code.start_codons, False

    starts = frozenset(normalise(codon) for codon in spec)
    if not starts:
        raise ValueError("at least one start codon is required")
    for codon in sorted(starts):
        if len(codon) != 3 or set(codon) - set("ACGT"):
            raise ValueError(f"not a codon: {codon!r}")
    return starts, False


def find_orfs(
    sequence: str,
    min_length: int = 300,
    table: int | str | GeneticCode = STANDARD_CODE,
    start_codons: Sequence[str] | Iterable[str] | str | None = "atg",
    strand: str = "both",
    include_truncated: bool = False,
    ignore_nested: bool = False,
    offset: int = 0,
) -> list[ORF]:
    """Find ORFs in all requested reading frames.

    Args:
        sequence: nucleotide sequence; U is read as T and case is ignored.
        min_length: smallest ORF to report, in nucleotides including the stop codon.
        table: NCBI translation table number, name, alias or ``GeneticCode``.
        start_codons: one of the modes in :data:`START_MODES`, or an explicit
            iterable of codons.
        strand: 'both', 'plus' or 'minus'.
        include_truncated: also report an ORF that runs off the end of the
            sequence without reaching a stop codon.
        ignore_nested: drop every ORF completely placed within a longer ORF on
            the same strand, applied after ``min_length``.
        offset: number of bases preceding ``sequence`` in the enclosing
            sequence, added to every reported coordinate.

    Returns:
        ORFs sorted by descending length, then start coordinate, then frame,
        numbered in that order.

    Raises:
        ValueError: on an unknown strand or start mode, a negative
            ``min_length``, or a malformed explicit start codon.
        KeyError: on an unknown translation table.
    """
    if strand not in STRANDS:
        raise ValueError(f"strand must be one of {', '.join(STRANDS)}; got {strand!r}")
    if min_length < 0:
        raise ValueError(f"min_length must not be negative; got {min_length}")

    code = table if isinstance(table, GeneticCode) else get_genetic_code(table)
    starts, any_sense = _resolve_starts(start_codons, code)
    stops = code.stop_codons

    work = normalise(sequence)
    total = len(work)

    hits: list[tuple[int, tuple[int, int, str, str | None, str, bool, bool]]] = []
    if strand in ("both", "plus"):
        hits.extend(
            _scan_strand(work, total, 1, starts, stops, any_sense, include_truncated)
        )
    if strand in ("both", "minus"):
        hits.extend(
            _scan_strand(
                reverse_complement(work), total, -1, starts, stops, any_sense,
                include_truncated,
            )
        )

    kept = [hit for hit in hits if len(hit[1][4]) >= min_length]
    kept.sort(key=lambda hit: (-len(hit[1][4]), hit[1][0], hit[0]))

    orfs: list[ORF] = []
    for frame, (start, end, start_codon, stop_codon, dna, truncated, partial) in kept:
        has_initiation_codon = not any_sense and start_codon in code.start_codons
        peptide = translate(dna, code)
        if stop_codon is not None:
            peptide = peptide[:-1]
        # An initiation codon is translated as methionine, which is how every
        # NCBI table marks its starts and how deposited CDS translations read.
        if has_initiation_codon and peptide:
            peptide = "M" + peptide[1:]
        orfs.append(
            ORF(
                index=0,
                strand="+" if frame > 0 else "-",
                frame=frame,
                start=start + offset,
                end=end + offset,
                length_nt=len(dna),
                start_codon=start_codon,
                stop_codon=stop_codon,
                protein=peptide,
                truncated=truncated,
                five_prime_partial=partial,
                has_initiation_codon=has_initiation_codon,
            )
        )

    if ignore_nested:
        orfs = _drop_nested(orfs)

    return [replace(orf, index=number) for number, orf in enumerate(orfs, start=1)]
