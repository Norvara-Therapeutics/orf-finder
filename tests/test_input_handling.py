"""Input quirks and output-path safety."""

import io

import pytest

from orf_finder.cli import _unique_names, main
from orf_finder.fasta import FastaError, Record, parse_fasta, read_fasta

BOM = "\ufeff"


def _run(argv):
    out = io.StringIO()
    code = main(argv, stdout=out)
    return code, out.getvalue()


def test_a_byte_order_mark_does_not_hide_the_first_header():
    (record,) = list(parse_fasta(io.StringIO(BOM + ">s desc\nACGT\n")))
    assert record.identifier == "s"
    assert record.description == "desc"
    assert record.sequence == "ACGT"


def test_a_byte_order_mark_is_stripped_when_reading_a_file(tmp_path):
    path = tmp_path / "bom.fasta"
    path.write_text(BOM + ">s\nATGAAACCCGGGTAA\n", encoding="utf-8")
    (record,) = read_fasta(path)
    assert record.identifier == "s"
    code, out = _run(["scan", "-i", str(path), "-m", "15", "-f", "tsv"])
    assert code == 0
    assert "MKPG" in out


def test_a_byte_order_mark_only_counts_on_the_first_line():
    # Mid-file it is not a mark, so it must not be silently absorbed.
    with pytest.raises(FastaError, match="not nucleotide sequence"):
        list(parse_fasta(io.StringIO(">s\nACGT\n" + BOM + ">t\nACGT\n")))


def test_header_tabs_separate_identifier_from_description(tmp_path):
    path = tmp_path / "tab.fasta"
    path.write_text(">seq\tDescription text\nATGAAATAA\n")
    (record,) = read_fasta(path)
    assert (record.identifier, record.description) == ("seq", "Description text")
    _, output = _run(["scan", "-i", str(path), "-m", "9", "--strand", "plus", "-f", "tsv"])
    header, row = output.splitlines()
    assert len(row.split("\t")) == len(header.split("\t"))


@pytest.mark.parametrize(
    "text",
    [
        ">s\nHello world, this is not DNA\n",
        ">s\nATGAAA>oops\nCCCGGGTAA\n",
        ">s\n1 2 3 4 5\n",
        ">s\nMKPGSTOP*\n",
    ],
)
def test_text_that_is_not_a_sequence_is_refused(text):
    """Prose or a stray marker must not be measured as though it were DNA."""
    with pytest.raises(FastaError, match="not nucleotide sequence"):
        list(parse_fasta(io.StringIO(text)))


def test_the_error_names_the_line_and_the_offending_characters():
    with pytest.raises(FastaError) as error:
        list(parse_fasta(io.StringIO(">s\nACGT\nACGZ\n")))
    assert "line 3" in str(error.value)
    assert "'Z'" in str(error.value)


@pytest.mark.parametrize(
    "sequence",
    [
        "ACGT",          # unambiguous DNA
        "acgt",          # soft-masked
        "ACGU",          # RNA
        "RYSWKMBDHVN",   # IUPAC ambiguity codes
        "ATG---AAATAA",  # alignment gaps
        "ATG...AAATAA",
    ],
)
def test_valid_nucleotide_text_is_accepted(sequence):
    (record,) = list(parse_fasta(io.StringIO(f">s\n{sequence}\n")))
    assert record.sequence == sequence


def test_a_directory_is_reported_as_such(tmp_path):
    with pytest.raises(SystemExit) as exit_info:
        _run(["scan", "-i", str(tmp_path)])
    assert exit_info.value.code == 2


@pytest.mark.parametrize(
    ("identifiers", "expected"),
    [
        (["a", "b"], ["a", "b"]),
        (["a", "a"], ["a", "a.2"]),
        (["a", "a", "a"], ["a", "a.2", "a.3"]),
        (["a", "b", "a"], ["a", "b", "a.2"]),
        # A generated name must not collide with one already present.
        (["a", "a.2", "a"], ["a", "a.2", "a.3"]),
    ],
)
def test_duplicate_identifiers_are_disambiguated(identifiers, expected):
    records = [Record(name, "", "ACGT") for name in identifiers]
    names = _unique_names(records)
    assert names == expected
    assert len(set(names)) == len(names)


def test_duplicate_identifiers_give_unique_protein_fasta_names(tmp_path):
    path = tmp_path / "dup.fasta"
    path.write_text(">dup\nCCATGAAAAAAAAAAAATAAGG\n>dup\nCCATGTTTTTTTTTTTTTAAGG\n")
    peptides = tmp_path / "d.faa"
    _run(["scan", "-i", str(path), "-m", "18", "-f", "none",
          "--protein-fasta", str(peptides)])
    deflines = [line for line in peptides.read_text().split("\n") if line.startswith(">")]
    assert len(deflines) == 2
    assert len({line.split()[0] for line in deflines}) == 2


def test_feature_table_refuses_duplicate_source_identifiers(tmp_path):
    path = tmp_path / "dup.fasta"
    path.write_text(">dup\nATGAAATAA\n>dup\nATGCCCTAA\n")
    table = tmp_path / "d.tbl"
    with pytest.raises(SystemExit) as exit_info:
        _run(["scan", "-i", str(path), "-m", "9", "-f", "none",
              "--feature-table", str(table)])
    assert exit_info.value.code == 2
    assert not table.exists()


def test_the_tsv_keeps_the_original_identifier(tmp_path):
    path = tmp_path / "dup.fasta"
    path.write_text(">dup\nCCATGAAAAAAAAAAAATAAGG\n>dup\nCCATGTTTTTTTTTTTTTAAGG\n")
    _, out = _run(["scan", "-i", str(path), "-m", "18", "-f", "tsv"])
    records = [line.split("\t")[0] for line in out.strip().split("\n")[1:]]
    assert records == ["dup", "dup"]


def test_two_outputs_may_not_share_a_path(tmp_path):
    path = tmp_path / "in.fasta"
    path.write_text(">s\nCCATGAAAAAAAAAAAATAAGG\n")
    shared = str(tmp_path / "both.out")
    with pytest.raises(SystemExit) as exit_info:
        _run(["scan", "-i", str(path), "-f", "none", "-o", shared,
              "--protein-fasta", shared])
    assert exit_info.value.code == 2


def test_two_outputs_may_not_share_a_hard_link(tmp_path):
    path = tmp_path / "in.fasta"
    path.write_text(">s\nATGAAATAA\n")
    table = tmp_path / "orfs.tsv"
    table.write_text("original")
    peptides = tmp_path / "orfs.faa"
    peptides.hardlink_to(table)
    with pytest.raises(SystemExit) as exit_info:
        _run(["scan", "-i", str(path), "-m", "9", "-f", "none",
              "-o", str(table), "--protein-fasta", str(peptides)])
    assert exit_info.value.code == 2
    assert table.read_text() == "original"


def test_output_may_not_overwrite_input_even_through_a_hard_link(tmp_path):
    path = tmp_path / "in.fasta"
    original = ">s\nATGAAATAA\n"
    path.write_text(original)
    alias = tmp_path / "alias.tsv"
    alias.hardlink_to(path)
    for output in (path, alias):
        with pytest.raises(SystemExit) as exit_info:
            _run(["scan", "-i", str(path), "-m", "9", "-f", "none", "-o", str(output)])
        assert exit_info.value.code == 2
        assert path.read_text() == original


def test_distinct_output_paths_are_accepted(tmp_path):
    path = tmp_path / "in.fasta"
    path.write_text(">s\nCCATGAAAAAAAAAAAATAAGG\n")
    code, _ = _run(["scan", "-i", str(path), "-m", "18", "-f", "none",
                    "-o", str(tmp_path / "a.tsv"),
                    "--protein-fasta", str(tmp_path / "b.faa")])
    assert code == 0
    assert (tmp_path / "a.tsv").exists() and (tmp_path / "b.faa").exists()
