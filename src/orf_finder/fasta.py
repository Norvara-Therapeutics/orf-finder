"""Minimal FASTA reader.

Accepts plain text from any path or open handle, tolerates blank lines, CRLF
line endings, a leading byte-order mark and internal whitespace, and preserves
record order. Sequence characters are checked against the IUPAC nucleotide
alphabet so that a file which is not a sequence is refused rather than measured.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import IO, Iterable, Iterator


class FastaError(ValueError):
    """Raised when input cannot be parsed as FASTA."""


#: Characters permitted in a sequence: the IUPAC nucleotide codes, including
#: the ambiguity codes, plus the gap characters used by alignment formats.
VALID_RESIDUES: frozenset[str] = frozenset("ACGTURYSWKMBDHVN" + "acgturyswkmbdhvn" + "-.~")


def _check_residues(line: str, lineno: int) -> None:
    """Reject sequence data that is not nucleotide text.

    Without this, a file that is not a sequence at all — prose, a saved error
    page, a column of numbers — is absorbed as residues and reported with a
    composition and a molecular weight, which is worse than an error.
    """
    invalid = sorted({character for character in line if character not in VALID_RESIDUES})
    if invalid:
        shown = ", ".join(repr(character) for character in invalid[:5])
        if len(invalid) > 5:
            shown += ", ..."
        raise FastaError(
            f"line {lineno}: not nucleotide sequence; unexpected character(s): {shown}"
        )


@dataclass(frozen=True)
class Record:
    """One FASTA record.

    Attributes:
        identifier: text up to the first whitespace of the header line.
        description: remainder of the header line, empty when absent.
        sequence: residues with all whitespace removed, case preserved.
    """

    identifier: str
    description: str
    sequence: str

    @property
    def header(self) -> str:
        """The header line without its '>' prefix."""
        return f"{self.identifier} {self.description}".strip()


def parse_fasta(lines: Iterable[str]) -> Iterator[Record]:
    """Yield records from an iterable of lines."""
    identifier: str | None = None
    description = ""
    chunks: list[str] = []
    seen_header = False

    for lineno, raw in enumerate(lines, start=1):
        if lineno == 1:
            # A file exported on Windows may open with a byte-order mark, which
            # would otherwise make the first header unrecognisable.
            raw = raw.lstrip("\ufeff")
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if identifier is not None:
                yield Record(identifier, description, "".join(chunks))
            seen_header = True
            head = line[1:].strip()
            if not head:
                raise FastaError(f"line {lineno}: header has no identifier")
            fields = head.split(None, 1)
            identifier = fields[0]
            description = fields[1] if len(fields) > 1 else ""
            chunks = []
        else:
            if not seen_header:
                raise FastaError(f"line {lineno}: sequence data before any '>' header")
            residues = "".join(line.split())
            _check_residues(residues, lineno)
            chunks.append(residues)

    if identifier is not None:
        yield Record(identifier, description, "".join(chunks))
    elif not seen_header:
        raise FastaError("input contains no FASTA records")


def read_fasta(source: str | Path | IO[str]) -> list[Record]:
    """Read every record from a path or an open text handle."""
    if hasattr(source, "read"):
        return list(parse_fasta(source))  # type: ignore[arg-type]
    path = Path(source)
    with path.open("r", encoding="utf-8-sig") as handle:
        return list(parse_fasta(handle))
