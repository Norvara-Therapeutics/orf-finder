"""Command-line interface."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import IO, Sequence

from orf_finder import __version__
from orf_finder.fasta import FastaError, Record, read_fasta
from orf_finder.genetic_code import (
    available_tables,
    get_genetic_code,
    resolve_table_name,
)
from orf_finder.orf import ORF, START_MODES, STRANDS, cds_sequence, find_orfs
from orf_finder.report import RecordReport, render_html
from orf_finder.sequence import SequenceMetrics, metrics

STDOUT_FORMATS = ("table", "tsv", "json", "none")

_TSV_COLUMNS = (
    "record",
    "orf",
    "strand",
    "frame",
    "start",
    "end",
    "length_nt",
    "length_aa",
    "start_codon",
    "stop_codon",
    "truncated",
    "five_prime_partial",
    "has_initiation_codon",
    "protein",
)


def _start_codons(text: str) -> str | tuple[str, ...]:
    """Resolve --start-codons: a mode keyword or a comma-separated list of codons."""
    keyword = text.strip().lower()
    if keyword in START_MODES:
        return keyword
    if keyword in ("default", "alternative"):
        return "atg" if keyword == "default" else "alt"
    codons = tuple(
        codon.strip().upper().replace("U", "T")
        for codon in text.split(",")
        if codon.strip()
    )
    if not codons:
        raise argparse.ArgumentTypeError("no start codons given")
    for codon in codons:
        if len(codon) != 3 or set(codon) - set("ACGT"):
            raise argparse.ArgumentTypeError(f"not a codon: {codon!r}")
    return codons


def _table(text: str) -> int:
    """Resolve --table: an NCBI table number, name or alias."""
    try:
        return resolve_table_name(text)
    except KeyError:
        raise argparse.ArgumentTypeError(
            f"unknown translation table: {text!r}; "
            "run 'orf-finder tables' for the available tables"
        ) from None


def _positive(text: str) -> int:
    value = _non_negative(text)
    if value < 1:
        raise argparse.ArgumentTypeError(f"must be 1 or greater: {value}")
    return value


def _non_negative(text: str) -> int:
    try:
        value = int(text)
    except ValueError:
        raise argparse.ArgumentTypeError(f"not an integer: {text!r}") from None
    if value < 0:
        raise argparse.ArgumentTypeError(f"must not be negative: {value}")
    return value


def build_parser() -> argparse.ArgumentParser:
    """Construct the argument parser."""
    parser = argparse.ArgumentParser(
        prog="orf-finder",
        description=(
            "Compute nucleic-acid metrics and scan all six reading frames for "
            "open reading frames."
        ),
    )
    parser.add_argument("--version", action="version", version=f"orf-finder {__version__}")
    subparsers = parser.add_subparsers(dest="command", required=True)

    scan = subparsers.add_parser(
        "scan",
        help="scan a FASTA file for ORFs",
        description=(
            "Scan every record of a FASTA file. Coordinates are 1-based and "
            "inclusive on the forward strand, with start <= end on both strands; "
            "a minus-strand ORF is read from end towards start."
        ),
    )
    scan.add_argument(
        "-i",
        "--input",
        required=True,
        metavar="PATH",
        help="FASTA file to scan, or '-' to read standard input",
    )
    scan.add_argument(
        "-m",
        "--min-length",
        type=_non_negative,
        default=300,
        metavar="N",
        help=(
            "smallest ORF to report, in nucleotides including the stop codon "
            "(default: 300)"
        ),
    )
    scan.add_argument(
        "-g",
        "--table",
        type=_table,
        default=1,
        metavar="TABLE",
        help=(
            "NCBI translation table as a number, name or alias, e.g. 11, "
            "bacterial, vertebrate-mito (default: 1, the standard code)"
        ),
    )
    scan.add_argument(
        "-b",
        "--begin",
        type=_positive,
        default=1,
        metavar="N",
        help="first base of the sequence fragment to scan, 1-based (default: 1)",
    )
    scan.add_argument(
        "-e",
        "--end",
        type=_non_negative,
        default=0,
        metavar="N",
        help="last base of the fragment to scan; 0 means to the end (default: 0)",
    )
    scan.add_argument(
        "--strand",
        choices=STRANDS,
        default="both",
        help="strand(s) to scan (default: both)",
    )
    scan.add_argument(
        "--start-codons",
        type=_start_codons,
        default="atg",
        metavar="CODONS",
        help=(
            "start codons to use: 'atg' (default), 'alt' for ATG plus the "
            "translation table's alternative initiation codons, 'any' for any "
            "sense codon (every region between stops), or an explicit "
            "comma-separated list of codons"
        ),
    )
    scan.add_argument(
        "-n",
        "--ignore-nested",
        action="store_true",
        help="drop every ORF completely placed within a longer ORF on the same strand",
    )
    scan.add_argument(
        "--include-truncated",
        action="store_true",
        help="also report an ORF that reaches the end of the sequence with no stop codon",
    )
    scan.add_argument(
        "--report",
        metavar="PATH",
        help="write a self-contained HTML report to PATH",
    )
    scan.add_argument(
        "--no-peptides",
        action="store_true",
        help="omit translated peptide sequences from the HTML report",
    )
    scan.add_argument(
        "-o",
        "--output",
        metavar="PATH",
        help="write the ORF table to PATH as TSV",
    )
    scan.add_argument(
        "--protein-fasta",
        metavar="PATH",
        help="write translated ORF peptides to PATH as FASTA",
    )
    scan.add_argument(
        "--cds-fasta",
        metavar="PATH",
        help="write ORF coding sequences to PATH as nucleotide FASTA",
    )
    scan.add_argument(
        "--feature-table",
        metavar="PATH",
        help="write a five-column NCBI feature table to PATH (requires unique FASTA IDs)",
    )
    scan.add_argument(
        "-f",
        "--format",
        choices=STDOUT_FORMATS,
        default="table",
        dest="stdout_format",
        help="what to print on standard output (default: table)",
    )
    subparsers.add_parser(
        "tables",
        help="list the available NCBI translation tables",
        description="List every NCBI translation table bundled with orf-finder.",
    )
    return parser


def _tsv_rows(reports: Sequence[RecordReport]) -> list[list[str]]:
    rows = [list(_TSV_COLUMNS)]
    for item in reports:
        for orf in item.orfs:
            rows.append(
                [
                    item.record.identifier,
                    orf.label,
                    orf.strand,
                    str(orf.frame),
                    str(orf.start),
                    str(orf.end),
                    str(orf.length_nt),
                    str(orf.length_aa),
                    orf.start_codon,
                    orf.stop_codon or "",
                    "yes" if orf.truncated else "no",
                    "yes" if orf.five_prime_partial else "no",
                    "yes" if orf.has_initiation_codon else "no",
                    orf.protein,
                ]
            )
    return rows


def format_tsv(reports: Sequence[RecordReport]) -> str:
    """Render the ORF table as TSV, header included."""
    return "".join("\t".join(row) + "\n" for row in _tsv_rows(reports))


def _defline(identifier: str, orf: ORF, measure: str) -> str:
    """Build a FASTA defline shared by the peptide and coding-sequence outputs."""
    sign = "+" if orf.frame > 0 else ""
    partial_ends: list[str] = []
    if orf.five_prime_partial or not orf.has_initiation_codon:
        partial_ends.append("5prime")
    if orf.truncated:
        partial_ends.append("3prime")
    flags = f" partial={','.join(partial_ends)}" if partial_ends else ""
    return (
        f">{identifier}_{orf.label} frame={sign}{orf.frame} "
        f"location={orf.start}-{orf.end} strand={orf.strand} {measure}{flags}\n"
    )


def format_protein_fasta(reports: Sequence[RecordReport], width: int = 60) -> str:
    """Render every reported peptide as FASTA."""
    parts: list[str] = []
    for item in reports:
        for orf in item.orfs:
            parts.append(_defline(item.name, orf, f"length_aa={orf.length_aa}"))
            peptide = orf.protein
            for i in range(0, len(peptide), width):
                parts.append(peptide[i : i + width] + "\n")
    return "".join(parts)


def format_cds_fasta(reports: Sequence[RecordReport], width: int = 60) -> str:
    """Render the coding sequence of every reported ORF as nucleotide FASTA."""
    parts: list[str] = []
    for item in reports:
        for orf in item.orfs:
            parts.append(_defline(item.name, orf, f"length_nt={orf.length_nt}"))
            dna = cds_sequence(orf, item.record.sequence)
            for i in range(0, len(dna), width):
                parts.append(dna[i : i + width] + "\n")
    return "".join(parts)


def format_feature_table(reports: Sequence[RecordReport], table_id: int = 1) -> str:
    """Render the ORFs as a five-column feature table.

    Coordinates run 5' to 3', so a minus-strand feature descends. A partial end
    carries the usual '<' or '>' marker. A record with no reported ORFs is
    omitted rather than given an empty block. Nonstandard genetic codes are
    recorded with a transl_table qualifier on every CDS.
    """
    parts: list[str] = []
    for item in reports:
        if not item.orfs:
            continue
        parts.append(f">Feature {item.name}\n")
        for orf in item.orfs:
            if orf.strand == "+":
                first, last = orf.start, orf.end
            else:
                first, last = orf.end, orf.start
            head = (
                f"<{first}"
                if orf.five_prime_partial or not orf.has_initiation_codon
                else str(first)
            )
            tail = f">{last}" if orf.truncated else str(last)
            parts.append(f"{head}\t{tail}\tCDS\n")
            parts.append("\t\t\tproduct\thypothetical protein\n")
            if table_id != 1:
                parts.append(f"\t\t\ttransl_table\t{table_id}\n")
    return "".join(parts)


def format_json(
    reports: Sequence[RecordReport], parameters: dict[str, str]
) -> str:
    """Render metrics and ORFs as JSON."""
    payload = {
        "tool": "orf-finder",
        "version": __version__,
        "parameters": parameters,
        "records": [
            {
                "identifier": item.record.identifier,
                "description": item.record.description,
                "metrics": {
                    "length": item.metrics.length,
                    "molecule": item.metrics.molecule,
                    "counts": item.metrics.counts,
                    "ambiguous": item.metrics.ambiguous,
                    "gaps": item.metrics.gaps,
                    "gc_percent": item.metrics.gc_percent,
                    "molecular_weight": item.metrics.molecular_weight,
                    "molecular_weight_duplex": item.metrics.molecular_weight_duplex,
                },
                "orfs": [orf.as_dict() for orf in item.orfs],
            }
            for item in reports
        ],
    }
    return json.dumps(payload, indent=2) + "\n"


def _format_metrics_line(summary: SequenceMetrics, fragment: bool = False) -> str:
    counts = " ".join(
        f"{base}={summary.counts[base]:,}" for base in summary.base_order
    )
    gc = "n/a" if summary.gc_percent is None else f"{summary.gc_percent:.2f}%"
    mass = (
        "n/a"
        if summary.molecular_weight is None
        else f"{summary.molecular_weight:,.2f} g/mol"
    )
    scope = " (fragment)" if fragment else ""
    gaps = f" gaps={summary.gaps:,}" if summary.gaps else ""
    return (
        f"  {summary.length:,} nt {summary.molecule}{scope}  {counts}"
        f" N={summary.counts['N']:,}{gaps}  GC={gc}  mass={mass}"
    )


def format_table(reports: Sequence[RecordReport], fragment: bool = False) -> str:
    """Render a readable per-record summary for the terminal.

    Args:
        reports: one entry per scanned record.
        fragment: whether the metrics describe a fragment of each record rather
            than the whole of it.
    """
    lines: list[str] = []
    for item in reports:
        header = item.record.identifier
        if item.record.description:
            header = f"{header}  {item.record.description}"
        lines.append(f"> {header}")
        lines.append(_format_metrics_line(item.metrics, fragment))
        if not item.orfs:
            lines.append("  no ORFs met the reporting threshold")
            lines.append("")
            continue
        lines.append(
            f"  {len(item.orfs)} ORF{'' if len(item.orfs) == 1 else 's'}"
        )
        lines.append(
            f"    {'ORF':<8}{'FRAME':>6}{'START':>12}{'END':>12}"
            f"{'NT':>9}{'AA':>8}  {'START':<6}{'STOP':<6}"
        )
        for orf in item.orfs:
            sign = "+" if orf.frame > 0 else ""
            lines.append(
                f"    {orf.label:<8}{sign + str(orf.frame):>6}{orf.start:>12,}"
                f"{orf.end:>12,}{orf.length_nt:>9,}{orf.length_aa:>8,}  "
                f"{orf.start_codon:<6}{orf.stop_codon or '-':<6}"
            )
        lines.append("")
    return "\n".join(lines).rstrip("\n") + "\n" if lines else ""


def _read_records(source: str) -> list[Record]:
    if source == "-":
        return read_fasta(sys.stdin)
    path = Path(source)
    if path.is_dir():
        raise IsADirectoryError(f"not a file: {source}")
    if not path.is_file():
        raise FileNotFoundError(f"no such file: {source}")
    return read_fasta(path)


def _describe_starts(start_codons: str | tuple[str, ...]) -> str:
    """Describe the start-codon setting for the parameter summary."""
    if isinstance(start_codons, tuple):
        return ", ".join(start_codons)
    return {
        "atg": "ATG only",
        "alt": "ATG and the table's alternative initiation codons",
        "any": "any sense codon (regions between stops)",
    }[start_codons]


def _parameters(args: argparse.Namespace, code_label: str) -> dict[str, str]:
    fragment = (
        "whole sequence"
        if args.begin == 1 and args.end == 0
        else f"{args.begin:,} to {'end' if args.end == 0 else format(args.end, ',')}"
    )
    return {
        "Minimum length": f"{args.min_length:,} nt (including stop codon)",
        "Strand": args.strand,
        "Start codons": _describe_starts(args.start_codons),
        "Genetic code": code_label,
        "Nested ORFs": "dropped" if args.ignore_nested else "reported",
        "Truncated ORFs": "reported" if args.include_truncated else "excluded",
        "Fragment": fragment,
    }


def format_tables() -> str:
    """Render the bundled translation tables."""
    lines = ["TABLE  NAME"]
    for code in available_tables():
        lines.append(f"{code.id:>5}  {code.name}")
    return "\n".join(lines) + "\n"


def _fragment(sequence: str, begin: int, end: int) -> str:
    """Slice the requested subrange.

    A fragment that starts beyond the end of the sequence is empty rather than
    an error, which is also how an empty FASTA record is handled.
    """
    if end and begin > end:
        raise ValueError(
            f"interval start is expected to be smaller than stop; got [{begin}, {end}]"
        )
    stop = len(sequence) if end == 0 else min(end, len(sequence))
    return sequence[begin - 1 : stop] if begin <= stop else ""


def _write(path: str, text: str) -> None:
    Path(path).write_text(text, encoding="utf-8")


def _unique_names(records: Sequence[Record]) -> list[str]:
    """Name each record uniquely, suffixing only identifiers that repeat.

    FASTA and feature-table consumers generally require distinct sequence
    names, while a FASTA file may legitimately repeat one. The second and later
    occurrences become ``<identifier>.2``, ``<identifier>.3`` and so on; a
    record whose identifier appears once keeps it unchanged.
    """
    totals: dict[str, int] = {}
    for record in records:
        totals[record.identifier] = totals.get(record.identifier, 0) + 1
    taken = set(totals)
    seen: dict[str, int] = {}
    names: list[str] = []
    for record in records:
        identifier = record.identifier
        seen[identifier] = seen.get(identifier, 0) + 1
        if totals[identifier] == 1 or seen[identifier] == 1:
            names.append(identifier)
            continue
        suffix = seen[identifier]
        candidate = f"{identifier}.{suffix}"
        # A generated name must not collide with a record named that already.
        while candidate in taken:
            suffix += 1
            candidate = f"{identifier}.{suffix}"
        taken.add(candidate)
        names.append(candidate)
    return names


def _check_output_paths(args: argparse.Namespace) -> None:
    """Refuse to overwrite the input or write different outputs to one file.

    Each output is written independently, so sharing a path would silently
    leave only whichever was written last. A path to the input would destroy
    the source FASTA after it has been read.
    """
    requested = {
        "--report": args.report,
        "--output": args.output,
        "--protein-fasta": args.protein_fasta,
        "--cds-fasta": args.cds_fasta,
        "--feature-table": args.feature_table,
    }
    claimed: dict[Path, str] = {}
    input_path = None if args.input == "-" else Path(args.input).expanduser().resolve()
    for flag, value in requested.items():
        if value is None:
            continue
        resolved = Path(value).expanduser().resolve()
        if input_path is not None and (
            resolved == input_path
            or (resolved.exists() and input_path.exists() and resolved.samefile(input_path))
        ):
            raise ValueError(f"{flag} would overwrite the input FASTA file: {value}")
        for other_path, other_flag in claimed.items():
            if resolved == other_path or (
                resolved.exists() and other_path.exists() and resolved.samefile(other_path)
            ):
                raise ValueError(
                    f"{other_flag} and {flag} would both write to {value}"
                )
        claimed[resolved] = flag


def run_scan(args: argparse.Namespace, stdout: IO[str]) -> int:
    """Execute the scan subcommand. Returns a process exit code."""
    _check_output_paths(args)
    records = _read_records(args.input)
    code = get_genetic_code(args.table)
    names = _unique_names(records)
    if args.feature_table and len(set(record.identifier for record in records)) != len(records):
        raise ValueError(
            "--feature-table requires unique input FASTA identifiers so feature "
            "headers match the source sequences"
        )
    reports = []
    for record, name in zip(records, names):
        fragment = _fragment(record.sequence, args.begin, args.end)
        reports.append(
            RecordReport(
                record=record,
                display_id=name,
                metrics=metrics(fragment),
                orfs=find_orfs(
                    fragment,
                    min_length=args.min_length,
                    table=code,
                    start_codons=args.start_codons,
                    strand=args.strand,
                    include_truncated=args.include_truncated,
                    ignore_nested=args.ignore_nested,
                    offset=args.begin - 1,
                ),
            )
        )
    parameters = _parameters(args, f"{code.id} ({code.name})")

    if args.report:
        _write(
            args.report,
            render_html(
                reports,
                source=("standard input" if args.input == "-" else Path(args.input).name),
                parameters=parameters,
                version=__version__,
                show_peptides=not args.no_peptides,
            ),
        )
    if args.output:
        _write(args.output, format_tsv(reports))
    if args.protein_fasta:
        _write(args.protein_fasta, format_protein_fasta(reports))
    if args.cds_fasta:
        _write(args.cds_fasta, format_cds_fasta(reports))
    if args.feature_table:
        _write(args.feature_table, format_feature_table(reports, table_id=code.id))

    if args.stdout_format == "table":
        stdout.write(format_table(reports, fragment=args.begin > 1 or args.end > 0))
    elif args.stdout_format == "tsv":
        stdout.write(format_tsv(reports))
    elif args.stdout_format == "json":
        stdout.write(format_json(reports, parameters))
    return 0


def main(argv: Sequence[str] | None = None, stdout: IO[str] | None = None) -> int:
    """Entry point. Returns a process exit code."""
    parser = build_parser()
    args = parser.parse_args(argv)
    out = sys.stdout if stdout is None else stdout
    try:
        if args.command == "tables":
            out.write(format_tables())
            return 0
        if args.command == "scan":
            return run_scan(args, out)
    except (FastaError, FileNotFoundError, KeyError, OSError, ValueError) as error:
        parser.exit(2, f"orf-finder: error: {error}\n")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
