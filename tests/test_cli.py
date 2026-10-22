import io
import json

import pytest

from orf_finder.cli import build_parser, main

FASTA = (
    ">seq1 a test record\n"
    "CCATGAAAAAAAAAAAATAAGG\n"
    ">seq2 no orfs here\n"
    "ACGTACGTACGT\n"
)


@pytest.fixture()
def fasta_file(tmp_path):
    path = tmp_path / "input.fasta"
    path.write_text(FASTA)
    return path


def _run(argv):
    out = io.StringIO()
    code = main(argv, stdout=out)
    return code, out.getvalue()


def test_table_output_summarises_every_record(fasta_file):
    code, out = _run(["scan", "-i", str(fasta_file), "-m", "18"])
    assert code == 0
    assert "> seq1  a test record" in out
    assert "22 nt DNA" in out
    assert "ORF1" in out
    assert "no ORFs met the reporting threshold" in out


def test_tsv_output_has_a_header_and_one_row_per_orf(fasta_file):
    _, out = _run(["scan", "-i", str(fasta_file), "-m", "18", "-f", "tsv"])
    lines = out.strip().split("\n")
    assert lines[0].split("\t")[:6] == [
        "record", "orf", "strand", "frame", "start", "end",
    ]
    assert len(lines) == 2
    fields = lines[1].split("\t")
    assert fields[0] == "seq1"
    assert fields[4:8] == ["3", "20", "18", "5"]
    assert fields[-1] == "MKKKK"


def test_json_output_carries_metrics_parameters_and_orfs(fasta_file):
    _, out = _run(["scan", "-i", str(fasta_file), "-m", "18", "-f", "json"])
    payload = json.loads(out)
    assert payload["tool"] == "orf-finder"
    assert payload["parameters"]["Strand"] == "both"
    assert [record["identifier"] for record in payload["records"]] == ["seq1", "seq2"]
    first = payload["records"][0]
    assert first["metrics"]["length"] == 22
    assert first["metrics"]["counts"]["A"] == 15
    assert first["orfs"][0]["protein"] == "MKKKK"
    assert payload["records"][1]["orfs"] == []


def test_format_none_prints_nothing(fasta_file):
    _, out = _run(["scan", "-i", str(fasta_file), "-f", "none"])
    assert out == ""


def test_report_tsv_and_protein_fasta_files_are_written(tmp_path, fasta_file):
    report = tmp_path / "report.html"
    table = tmp_path / "orfs.tsv"
    peptides = tmp_path / "orfs.faa"
    code, _ = _run([
        "scan", "-i", str(fasta_file), "-m", "18", "-f", "none",
        "--report", str(report), "-o", str(table), "--protein-fasta", str(peptides),
    ])
    assert code == 0
    html = report.read_text()
    assert html.startswith("<!DOCTYPE html>")
    assert "seq1" in html and "seq2" in html
    assert "input.fasta" in html
    assert table.read_text().count("\n") == 2
    assert peptides.read_text() == (
        ">seq1_ORF1 frame=+3 location=3-20 strand=+ length_aa=5\nMKKKK\n"
    )


def test_the_report_is_byte_identical_across_runs(tmp_path, fasta_file):
    first = tmp_path / "a.html"
    second = tmp_path / "b.html"
    for path in (first, second):
        _run(["scan", "-i", str(fasta_file), "-f", "none", "--report", str(path)])
    assert first.read_bytes() == second.read_bytes()


def test_reads_standard_input(monkeypatch):
    monkeypatch.setattr("sys.stdin", io.StringIO(FASTA))
    _, out = _run(["scan", "-i", "-", "-m", "18", "-f", "tsv"])
    assert "seq1" in out


def test_strand_and_start_codon_options(tmp_path):
    path = tmp_path / "gtg.fasta"
    path.write_text(">g\nGTGAAAAAAAAAAAATAA\n")
    # GTG is not an initiation codon of the standard code ...
    _, out = _run(["scan", "-i", str(path), "-m", "18", "--start-codons", "alt",
                   "--strand", "plus", "-f", "tsv"])
    assert len(out.strip().split("\n")) == 1
    # ... but it is one in the bacterial code, where it initiates as methionine.
    _, out = _run(["scan", "-i", str(path), "-m", "18", "--start-codons", "alt",
                   "--table", "bacterial", "--strand", "plus", "-f", "tsv"])
    assert "MKKKK" in out


def test_truncated_orfs_require_the_flag(tmp_path):
    path = tmp_path / "trunc.fasta"
    path.write_text(">t\nATGAAAAAAAAAAAAAAA\n")
    _, out = _run(["scan", "-i", str(path), "-m", "18", "-f", "tsv"])
    assert len(out.strip().split("\n")) == 1
    _, out = _run(["scan", "-i", str(path), "-m", "18", "-f", "tsv",
                   "--include-truncated"])
    assert "yes" in out.strip().split("\n")[1]


def test_tables_subcommand_lists_every_table():
    code, out = _run(["tables"])
    assert code == 0
    lines = out.strip().split("\n")
    assert lines[0].split() == ["TABLE", "NAME"]
    assert len(lines) == 28  # header plus 27 tables
    assert "    1  Standard" in out
    assert "Bacterial" in out and "Karyorelict Nuclear" in out


def test_table_option_changes_translation(tmp_path):
    path = tmp_path / "mito.fasta"
    # TGA is a stop under the standard code and tryptophan in vertebrate mitochondria.
    path.write_text(">m\nATGAAATGATTTCCCGGGTAA\n")
    _, standard = _run(["scan", "-i", str(path), "-m", "9", "-f", "tsv", "--strand", "plus"])
    assert "MK" in standard
    _, mito = _run(["scan", "-i", str(path), "-m", "21", "-f", "tsv",
                    "--table", "vertebrate-mito", "--strand", "plus"])
    assert "MKWFPG" in mito


def test_feature_table_records_nonstandard_translation_table(tmp_path):
    path = tmp_path / "mito.fasta"
    path.write_text(">m\nATGAAATGATTTCCCGGGTAA\n")
    table = tmp_path / "mito.tbl"
    _run(["scan", "-i", str(path), "-m", "21", "--strand", "plus",
          "--table", "2", "--feature-table", str(table), "-f", "none"])
    assert "\t\t\ttransl_table\t2\n" in table.read_text()


def test_unknown_table_exits_with_code_two(fasta_file):
    with pytest.raises(SystemExit) as exit_info:
        _run(["scan", "-i", str(fasta_file), "--table", "martian"])
    assert exit_info.value.code == 2


def test_any_sense_mode_marks_five_prime_partial_orfs(tmp_path):
    path = tmp_path / "any.fasta"
    path.write_text(">a\nAAACCCTAAATGAAACCCGGGTAA\n")
    _, out = _run(["scan", "-i", str(path), "-m", "6", "-f", "tsv",
                   "--start-codons", "any", "--strand", "plus"])
    rows = [line.split("\t") for line in out.strip().split("\n")[1:]]
    header = _run(["scan", "-i", str(path), "-m", "6", "-f", "tsv"])[1].split("\n")[0]
    column = header.split("\t").index("five_prime_partial")
    assert any(row[column] == "yes" for row in rows)


def test_ignore_nested_option(tmp_path):
    path = tmp_path / "nested.fasta"
    outer = "ATG" + "AAA" * 40 + "TAA"
    sequence = outer[:4] + "ATG" + "CCC" * 10 + "TGA" + outer[4 + 36:]
    path.write_text(f">n\n{sequence}\n")
    _, before = _run(["scan", "-i", str(path), "-m", "33", "-f", "tsv", "--strand", "plus"])
    _, after = _run(["scan", "-i", str(path), "-m", "33", "-f", "tsv", "--strand", "plus",
                     "--ignore-nested"])
    assert len(before.strip().split("\n")) > len(after.strip().split("\n"))


def test_begin_and_end_restrict_the_scan_and_keep_absolute_coordinates(tmp_path):
    path = tmp_path / "frag.fasta"
    path.write_text(">f\n" + "TTTTT" + "ATGAAACCCGGGTAA" + "TTTTT" + "\n")
    _, whole = _run(["scan", "-i", str(path), "-m", "15", "-f", "tsv", "--strand", "plus"])
    assert whole.strip().split("\n")[1].split("\t")[4:6] == ["6", "20"]
    _, fragment = _run(["scan", "-i", str(path), "-m", "15", "-f", "tsv",
                        "--strand", "plus", "-b", "6", "-e", "20"])
    assert fragment.strip().split("\n")[1].split("\t")[4:6] == ["6", "20"]
    _, outside = _run(["scan", "-i", str(path), "-m", "15", "-f", "tsv",
                       "--strand", "plus", "-b", "7"])
    assert len(outside.strip().split("\n")) == 1


def test_an_empty_record_and_an_out_of_range_fragment_are_not_errors(tmp_path):
    path = tmp_path / "empty.fasta"
    path.write_text(">blank\n>seq1\nCCATGAAAAAAAAAAAATAAGG\n")
    code, out = _run(["scan", "-i", str(path), "-m", "18"])
    assert code == 0
    assert "> blank" in out and "ORF1" in out
    code, out = _run(["scan", "-i", str(path), "-m", "18", "-b", "500", "-f", "tsv"])
    assert code == 0
    assert len(out.strip().split("\n")) == 1


def test_an_inverted_interval_exits_with_code_two(fasta_file):
    with pytest.raises(SystemExit) as exit_info:
        _run(["scan", "-i", str(fasta_file), "-b", "20", "-e", "5"])
    assert exit_info.value.code == 2


def test_cds_fasta_translates_to_the_reported_peptide(tmp_path, fasta_file):
    cds = tmp_path / "orfs.fna"
    _run(["scan", "-i", str(fasta_file), "-m", "18", "-f", "none", "--cds-fasta", str(cds)])
    text = cds.read_text()
    assert text.startswith(">seq1_ORF1 frame=+3 location=3-20 strand=+ length_nt=18")
    assert "ATGAAAAAAAAAAAATAA" in text


def test_feature_table_uses_five_columns_and_descends_on_the_minus_strand(tmp_path):
    path = tmp_path / "both.fasta"
    path.write_text(">p\nCCATGAAAAAAAAAAAATAAGG\n>m\nCCTTATTTTTTTTTTTTCATGG\n")
    table = tmp_path / "orfs.tbl"
    _run(["scan", "-i", str(path), "-m", "18", "-f", "none", "--feature-table", str(table)])
    lines = table.read_text().split("\n")
    assert lines[0] == ">Feature p"
    assert lines[1] == "3\t20\tCDS"
    assert lines[2] == "\t\t\tproduct\thypothetical protein"
    assert lines[3] == ">Feature m"
    assert lines[4].split("\t")[:2] == ["20", "3"], "minus-strand features descend"


def test_feature_table_omits_a_record_with_no_orfs(tmp_path):
    path = tmp_path / "mixed.fasta"
    path.write_text(">blank\nACGTACGT\n>has\nCCATGAAAAAAAAAAAATAAGG\n")
    table = tmp_path / "m.tbl"
    _run(["scan", "-i", str(path), "-m", "18", "-f", "none", "--feature-table", str(table)])
    assert ">Feature blank" not in table.read_text()
    assert table.read_text().startswith(">Feature has")


def test_a_fragment_scan_says_so_in_the_terminal_output(tmp_path):
    path = tmp_path / "frag.fasta"
    path.write_text(">s\nGGGGGGGGGATGAAACCCGGGTAATTTTT\n")
    _, whole = _run(["scan", "-i", str(path), "-m", "15"])
    assert "29 nt DNA " in whole and "(fragment)" not in whole
    _, part = _run(["scan", "-i", str(path), "-m", "15", "-b", "10"])
    assert "20 nt DNA (fragment)" in part


def test_feature_table_marks_a_partial_end(tmp_path):
    path = tmp_path / "trunc.fasta"
    path.write_text(">t\nATGAAAAAAAAAAAAAAA\n")
    table = tmp_path / "t.tbl"
    _run(["scan", "-i", str(path), "-m", "18", "-f", "none", "--include-truncated",
          "--feature-table", str(table), "--strand", "plus"])
    assert table.read_text().split("\n")[1] == "1\t>18\tCDS"


def test_any_sense_without_initiator_marks_five_prime_end_in_outputs(tmp_path):
    path = tmp_path / "any.fasta"
    path.write_text(">s\nTAACCCTAA\n")
    table = tmp_path / "any.tbl"
    cds = tmp_path / "any.fna"
    _, output = _run(["scan", "-i", str(path), "-m", "6", "--strand", "plus",
                      "--start-codons", "any", "-f", "tsv",
                      "--feature-table", str(table), "--cds-fasta", str(cds)])
    header, row = (line.split("\t") for line in output.splitlines())
    assert row[header.index("has_initiation_codon")] == "no"
    assert table.read_text().splitlines()[1] == "<4\t9\tCDS"
    assert "partial=5prime" in cds.read_text().splitlines()[0]


def test_doubly_partial_fasta_defline_has_one_partial_field(tmp_path):
    path = tmp_path / "partial.fasta"
    path.write_text(">p\nAAACCC\n")
    peptides = tmp_path / "partial.faa"
    _run(["scan", "-i", str(path), "-m", "6", "--strand", "plus",
          "--start-codons", "any", "--include-truncated", "-f", "none",
          "--protein-fasta", str(peptides)])
    defline = peptides.read_text().splitlines()[0]
    assert defline.count("partial=") == 1
    assert "partial=5prime,3prime" in defline


def test_missing_input_exits_with_code_two(tmp_path):
    with pytest.raises(SystemExit) as exit_info:
        _run(["scan", "-i", str(tmp_path / "absent.fasta")])
    assert exit_info.value.code == 2


def test_malformed_fasta_exits_with_code_two(tmp_path):
    path = tmp_path / "bad.fasta"
    path.write_text("ACGT\n")
    with pytest.raises(SystemExit) as exit_info:
        _run(["scan", "-i", str(path)])
    assert exit_info.value.code == 2


@pytest.mark.parametrize("argv", [
    ["scan", "-i", "x", "--start-codons", "AT"],
    ["scan", "-i", "x", "-m", "-5"],
    ["scan", "-i", "x", "--strand", "either"],
    ["scan"],
    [],
])
def test_usage_errors_exit_with_code_two(argv):
    with pytest.raises(SystemExit) as exit_info:
        build_parser().parse_args(argv)
    assert exit_info.value.code == 2


def test_version_is_reported():
    with pytest.raises(SystemExit) as exit_info:
        build_parser().parse_args(["--version"])
    assert exit_info.value.code == 0
