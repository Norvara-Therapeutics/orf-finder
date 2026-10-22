import io

import pytest

from orf_finder.fasta import FastaError, parse_fasta, read_fasta


def test_parses_multiple_records_and_splits_the_description():
    text = ">one first record\nACGT\nACGT\n>two\nTTTT\n"
    records = list(parse_fasta(io.StringIO(text)))
    assert [record.identifier for record in records] == ["one", "two"]
    assert records[0].description == "first record"
    assert records[0].sequence == "ACGTACGT"
    assert records[1].description == ""
    assert records[1].header == "two"


def test_tolerates_blank_lines_crlf_and_internal_whitespace():
    text = ">one\r\n\r\nAC GT\r\n\tAC\tGT \r\n"
    (record,) = list(parse_fasta(io.StringIO(text)))
    assert record.sequence == "ACGTACGT"


def test_preserves_case():
    (record,) = list(parse_fasta(io.StringIO(">one\nAcGt\n")))
    assert record.sequence == "AcGt"


def test_allows_an_empty_record():
    records = list(parse_fasta(io.StringIO(">empty\n>two\nAC\n")))
    assert records[0].sequence == ""
    assert records[1].sequence == "AC"


def test_rejects_sequence_before_a_header():
    with pytest.raises(FastaError, match="before any"):
        list(parse_fasta(io.StringIO("ACGT\n>one\nACGT\n")))


def test_rejects_an_empty_header():
    with pytest.raises(FastaError, match="no identifier"):
        list(parse_fasta(io.StringIO(">\nACGT\n")))


def test_rejects_input_without_records():
    with pytest.raises(FastaError, match="no FASTA records"):
        list(parse_fasta(io.StringIO("\n\n")))


def test_reads_from_a_path(tmp_path):
    path = tmp_path / "seq.fasta"
    path.write_text(">one desc\nACGTACGT\n")
    (record,) = read_fasta(path)
    assert record.sequence == "ACGTACGT"


def test_reads_from_a_handle(tmp_path):
    path = tmp_path / "seq.fasta"
    path.write_text(">one\nACGT\n")
    with path.open() as handle:
        (record,) = read_fasta(handle)
    assert record.identifier == "one"
