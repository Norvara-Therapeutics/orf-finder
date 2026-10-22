"""Self-contained HTML reporting.

The document embeds its own style rules, uses no JavaScript and loads nothing
from the network, so a report stays readable from a local file or an email
attachment years later. Output is byte-deterministic for a given input and set
of parameters: nothing in it depends on the clock or the filesystem beyond the
source name that is passed in.
"""

from __future__ import annotations

from dataclasses import dataclass
from html import escape
from typing import Sequence

from orf_finder.fasta import Record
from orf_finder.orf import ORF
from orf_finder.sequence import SequenceMetrics

_STYLE = """
html { -webkit-text-size-adjust: 100%; }
body { font: 16px/1.55 Georgia, 'Times New Roman', serif; color: #1a1a1a;
       background: #fff; margin: 0 auto; padding: 3rem 1.25rem 4rem;
       max-width: 46rem; }
h1 { font-size: 1.5rem; font-weight: normal; margin: 0 0 .25rem; }
h2, h3 { font-size: .8rem; font-variant: small-caps; letter-spacing: .08em;
     font-weight: normal; margin: 2.5rem 0 .6rem; border-bottom: 1px solid #1a1a1a;
     padding-bottom: .2rem; }
h3 { margin-top: 1.8rem; }
h2.record { font-size: 1.1rem; font-variant: normal; letter-spacing: normal;
     font-style: italic; border-bottom: none; padding-bottom: 0; margin: 2.8rem 0 .2rem; }
p { margin: .5rem 0; }
.lede { color: #555; font-size: .9rem; margin-bottom: .5rem; }
.note { color: #555; font-size: .85rem; }
code, .seq, td.num, th.num { font-family: ui-monospace, 'SF Mono', Menlo, Consolas, monospace; }
table { border-collapse: collapse; width: 100%; margin: .6rem 0 .2rem;
        font-size: .9rem; }
th { text-align: left; font-weight: normal; font-variant: small-caps;
     letter-spacing: .04em; border-bottom: 1px solid #1a1a1a; padding: .25rem .5rem .25rem 0; }
td { padding: .22rem .5rem .22rem 0; border-bottom: 1px solid #e6e6e6;
     vertical-align: baseline; }
td.num, th.num { text-align: right; padding-right: 1rem; font-size: .85rem; }
tr:last-child td { border-bottom: none; }
dl.metrics { display: grid; grid-template-columns: max-content 1fr; gap: .15rem 1.2rem;
             margin: .6rem 0; font-size: .9rem; }
dl.metrics dt { color: #555; }
dl.metrics dd { margin: 0; font-family: ui-monospace, 'SF Mono', Menlo, Consolas, monospace; }
.seq { font-size: .8rem; line-height: 1.45; white-space: pre-wrap; word-break: break-all;
       margin: .2rem 0 1rem; }
.bar { display: block; height: 2px; background: #1a1a1a; margin: .3rem 0 .1rem; }
.tag { font-variant: small-caps; letter-spacing: .04em; font-size: .8rem; color: #0f6f6a; }
.tag.warn { color: #8a5a00; }
footer { margin-top: 3rem; padding-top: .5rem; border-top: 1px solid #e6e6e6;
         color: #555; font-size: .8rem; }
@media print { body { padding: 0; max-width: none; } h2, h3 { page-break-after: avoid; } }
"""


@dataclass(frozen=True)
class RecordReport:
    """Everything the report needs about one scanned record.

    ``display_id`` disambiguates a record whose identifier is repeated in the
    input for peptide and coding-sequence FASTA outputs. Feature tables require
    unique input identifiers. It defaults to the record's own identifier.
    """

    record: Record
    metrics: SequenceMetrics
    orfs: Sequence[ORF]
    display_id: str | None = None

    @property
    def name(self) -> str:
        """The identifier to use where sequence names must be unique."""
        return self.display_id or self.record.identifier


def _thousands(value: int) -> str:
    return f"{value:,}"


def _wrap(sequence: str, width: int = 60) -> str:
    return "\n".join(sequence[i : i + width] for i in range(0, len(sequence), width))


def _metrics_block(metrics: SequenceMetrics) -> list[str]:
    rows: list[tuple[str, str]] = [
        ("Length", f"{_thousands(metrics.length)} nt"),
        ("Molecule", metrics.molecule),
    ]
    total = sum(metrics.counts[base] for base in metrics.base_order)
    for base in metrics.base_order:
        count = metrics.counts[base]
        share = f" ({100.0 * count / total:.1f}%)" if total else ""
        rows.append((base, f"{_thousands(count)}{share}"))
    rows.append(("N", _thousands(metrics.counts["N"])))
    if metrics.ambiguous != metrics.counts["N"]:
        rows.append(("Other ambiguous", _thousands(metrics.ambiguous - metrics.counts["N"])))
    if metrics.gaps:
        rows.append(("Alignment gaps", _thousands(metrics.gaps)))
    rows.append(
        ("GC content", "n/a" if metrics.gc_percent is None else f"{metrics.gc_percent:.2f}%")
    )
    if metrics.molecular_weight is not None:
        rows.append(("Mass, single strand", f"{metrics.molecular_weight:,.2f} g/mol"))
    if metrics.molecular_weight_duplex is not None:
        rows.append(("Mass, duplex", f"{metrics.molecular_weight_duplex:,.2f} g/mol"))

    parts = ['<dl class="metrics">']
    for name, value in rows:
        parts.append(f"<dt>{escape(name)}</dt><dd>{escape(value)}</dd>")
    parts.append("</dl>")
    if metrics.ambiguous:
        parts.append(
            '<p class="note">Mass is approximate: ambiguous positions are charged the '
            "mean residue mass of the four standard bases. GC content is computed over "
            "unambiguous bases only.</p>"
        )
    if metrics.gaps:
        parts.append(
            '<p class="note">Alignment gaps count toward coordinate length but are '
            "excluded from base composition and mass.</p>"
        )
    return parts


def _frame_block(orfs: Sequence[ORF], strand: str) -> list[str]:
    frames = [1, 2, 3, -1, -2, -3]
    if strand == "plus":
        frames = frames[:3]
    elif strand == "minus":
        frames = frames[3:]
    counts = {frame: 0 for frame in frames}
    for orf in orfs:
        counts[orf.frame] += 1
    widest = max(counts.values()) or 1

    parts = ["<table>", '<tr><th>Frame</th><th class="num">ORFs</th><th></th></tr>']
    for frame in frames:
        count = counts[frame]
        width = 0.0 if not count else max(1.0, 100.0 * count / widest)
        bar = f'<span class="bar" style="width: {width:.1f}%"></span>' if count else ""
        sign = "+" if frame > 0 else ""
        parts.append(
            f"<tr><td>{sign}{frame}</td>"
            f'<td class="num">{count}</td><td>{bar}</td></tr>'
        )
    parts.append("</table>")
    return parts


def _orf_table(orfs: Sequence[ORF]) -> list[str]:
    parts = [
        "<table>",
        "<tr><th>ORF</th><th>Strand</th><th>Frame</th>"
        '<th class="num">Start</th><th class="num">End</th>'
        '<th class="num">Length (nt)</th><th class="num">Length (aa)</th>'
        "<th>Start</th><th>Stop</th></tr>",
    ]
    for orf in orfs:
        stop = orf.stop_codon or '<span class="tag warn">none</span>'
        start = escape(orf.start_codon)
        if orf.five_prime_partial:
            start = f'{start} <span class="tag warn">5&prime;</span>'
        if not orf.has_initiation_codon:
            start = f'{start} <span class="tag warn">no initiator</span>'
        sign = "+" if orf.frame > 0 else ""
        parts.append(
            f"<tr><td>{escape(orf.label)}</td><td>{escape(orf.strand)}</td>"
            f"<td>{sign}{orf.frame}</td>"
            f'<td class="num">{_thousands(orf.start)}</td>'
            f'<td class="num">{_thousands(orf.end)}</td>'
            f'<td class="num">{_thousands(orf.length_nt)}</td>'
            f'<td class="num">{_thousands(orf.length_aa)}</td>'
            f"<td><code>{start}</code></td><td><code>{stop}</code></td></tr>"
        )
    parts.append("</table>")
    if any(orf.truncated for orf in orfs):
        parts.append(
            '<p class="note">An ORF with no stop codon runs off the end of the '
            "sequence; its length counts only complete codons.</p>"
        )
    if any(orf.five_prime_partial for orf in orfs):
        parts.append(
            '<p class="note">An ORF marked 5&prime; is the first region of its '
            "frame and has no upstream stop codon, so its 5&prime; end is a region "
            "boundary rather than an initiation codon.</p>"
        )
    if any(not orf.has_initiation_codon for orf in orfs):
        parts.append(
            '<p class="note">In any-sense mode, region boundaries are not '
            "initiation codons; these ORFs are incomplete at the 5&prime; end.</p>"
        )
    return parts


def _record_section(item: RecordReport, show_peptides: bool, strand: str) -> list[str]:
    record = item.record
    parts = [f'<h2 class="record">{escape(record.identifier)}</h2>']
    if record.description:
        parts.append(f'<p class="lede">{escape(record.description)}</p>')

    parts.append("<h3>Sequence metrics</h3>")
    parts.extend(_metrics_block(item.metrics))

    parts.append("<h3>Reading frames</h3>")
    if not item.orfs:
        parts.append(
            '<p>No open reading frame met the reporting threshold in the scanned frames.</p>'
        )
        return parts
    parts.extend(_frame_block(item.orfs, strand))

    parts.append("<h3>Open reading frames</h3>")
    parts.extend(_orf_table(item.orfs))

    if show_peptides:
        parts.append("<h3>Predicted peptides</h3>")
        for orf in item.orfs:
            sign = "+" if orf.frame > 0 else ""
            parts.append(
                f'<p class="note">{escape(orf.label)} &middot; frame {sign}{orf.frame} '
                f"&middot; {_thousands(orf.start)}&ndash;{_thousands(orf.end)} "
                f"&middot; {orf.length_aa} aa</p>"
            )
            parts.append(f'<div class="seq">{escape(_wrap(orf.protein))}</div>')
    return parts


def render_html(
    reports: Sequence[RecordReport],
    source: str,
    parameters: dict[str, str],
    version: str,
    show_peptides: bool = True,
) -> str:
    """Render a complete, self-contained HTML document.

    Args:
        reports: one entry per scanned record, in input order.
        source: name of the input shown in the header.
        parameters: scan settings, rendered as an ordered summary list.
        version: tool version recorded in the footer.
        show_peptides: include the translated peptide of every reported ORF.
    """
    total_orfs = sum(len(item.orfs) for item in reports)
    plural = "" if len(reports) == 1 else "s"

    parts = [
        "<!DOCTYPE html>",
        '<html lang="en">',
        "<head>",
        '<meta charset="utf-8">',
        '<meta name="viewport" content="width=device-width, initial-scale=1">',
        f"<title>ORF scan of {escape(source)}</title>",
        f"<style>{_STYLE}</style>",
        "</head>",
        "<body>",
        "<h1>Open reading frame scan</h1>",
        f'<p class="lede">{escape(source)} &middot; {len(reports)} record{plural} '
        f"&middot; {_thousands(total_orfs)} ORFs reported</p>",
        "<h2>Scan parameters</h2>",
        '<dl class="metrics">',
    ]
    for name, value in parameters.items():
        parts.append(f"<dt>{escape(name)}</dt><dd>{escape(value)}</dd>")
    parts.append("</dl>")

    for item in reports:
        parts.extend(_record_section(item, show_peptides, parameters.get("Strand", "both")))

    parts.extend(
        [
            "<footer>",
            "<p>Coordinates are 1-based and inclusive on the forward strand of the "
            "input, with start &le; end on both strands; a minus-strand ORF is read "
            "from end towards start. Nucleotide lengths include the stop codon, which "
            "the peptide excludes. The translation table and the start-codon rule are "
            "listed under Scan parameters; where a start codon is required, one ORF is "
            "reported per in-frame stop codon, beginning at the 5&prime;-most start "
            "codon after the preceding stop.</p>",
            f"<p>orf-finder {escape(version)}</p>",
            "</footer>",
            "</body>",
            "</html>",
            "",
        ]
    )
    return "\n".join(parts)
