"""Sequence metrics and six-frame ORF scanning with no third-party dependencies."""

from orf_finder.genetic_code import (
    GeneticCode,
    available_tables,
    get_genetic_code,
    resolve_table_name,
)
from orf_finder.orf import ORF, START_MODES, cds_sequence, find_orfs
from orf_finder.sequence import SequenceMetrics, metrics, reverse_complement
from orf_finder.fasta import Record, read_fasta

__version__ = "1.1.3rc3"

__all__ = [
    "GeneticCode",
    "ORF",
    "Record",
    "START_MODES",
    "SequenceMetrics",
    "__version__",
    "available_tables",
    "cds_sequence",
    "find_orfs",
    "get_genetic_code",
    "metrics",
    "read_fasta",
    "resolve_table_name",
    "reverse_complement",
]
