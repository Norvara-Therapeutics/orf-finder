from orf_finder.fasta import Record
from orf_finder.orf import find_orfs
from orf_finder.report import RecordReport, render_html
from orf_finder.sequence import metrics

SEQUENCE = "CC" + "ATG" + "AAA" * 4 + "TAA" + "GG"
PARAMETERS = {"Minimum length": "18 nt", "Strand": "both"}


def _reports(identifier="seq1", description="a test record", sequence=SEQUENCE):
    record = Record(identifier, description, sequence)
    return [
        RecordReport(record, metrics(sequence), find_orfs(sequence, min_length=18))
    ]


def _render(**kwargs):
    return render_html(
        _reports(**kwargs), source="seq.fasta", parameters=PARAMETERS, version="1.0.0"
    )


def test_document_is_well_formed_and_self_contained():
    html = _render()
    assert html.startswith("<!DOCTYPE html>")
    assert html.rstrip().endswith("</html>")
    assert "<script" not in html.lower()
    assert "http://" not in html and "https://" not in html
    assert "<link" not in html.lower()
    assert "<style>" in html


def test_document_reports_the_metrics_parameters_and_orfs():
    html = _render()
    assert "seq.fasta" in html
    assert "Minimum length" in html and "18 nt" in html
    assert "22 nt" in html  # sequence length
    assert "22.73%" in html  # GC content over unambiguous bases
    assert "ORF1" in html
    assert "MKKKK" in html
    assert ">3<" in html and ">20<" in html  # inclusive coordinates


def test_peptides_can_be_omitted():
    html = render_html(
        _reports(), source="s", parameters=PARAMETERS, version="1.0.0",
        show_peptides=False,
    )
    assert "Predicted peptides" not in html
    assert "MKKKK" not in html


def test_identifiers_and_descriptions_are_escaped():
    html = _render(identifier="<script>x</script>", description="a & b <b>")
    assert "<script>" not in html
    assert "&lt;script&gt;" in html
    assert "a &amp; b &lt;b&gt;" in html


def test_rendering_is_deterministic():
    assert _render() == _render()


def test_a_record_without_orfs_is_stated_explicitly():
    html = render_html(
        [RecordReport(Record("empty", "", "ACGT"), metrics("ACGT"), [])],
        source="s", parameters=PARAMETERS, version="1.0.0",
    )
    assert "No open reading frame met the reporting threshold" in html
    assert "Open reading frames" not in html


def test_report_only_counts_scanned_frames():
    sequence = "ATGAAATAA"
    report = [RecordReport(Record("one", "", sequence), metrics(sequence),
                           find_orfs(sequence, min_length=9, strand="plus"))]
    html = render_html(report, source="s", parameters={"Strand": "plus"}, version="1")
    assert "<td>+1</td>" in html
    assert "<td>-1</td>" not in html


def test_no_orf_report_does_not_claim_six_frames_for_plus_scan():
    html = render_html(
        [RecordReport(Record("empty", "", "ACGT"), metrics("ACGT"), [])],
        source="s", parameters={"Strand": "plus"}, version="1",
    )
    assert "scanned frames" in html
    assert "six frames" not in html


def test_the_footer_defers_to_the_listed_parameters():
    html = _render()
    assert "standard genetic code" not in html
    assert "listed under Scan parameters" in html


def test_a_five_prime_partial_orf_is_marked_and_explained():
    sequence = "AAACCCTAAATGAAACCCGGGTAA"
    record = Record("any", "", sequence)
    orfs = find_orfs(sequence, min_length=6, start_codons="any", strand="plus")
    assert any(orf.five_prime_partial for orf in orfs)
    html = render_html(
        [RecordReport(record, metrics(sequence), orfs)],
        source="s", parameters=PARAMETERS, version="1.1.0",
    )
    assert "5&prime;</span>" in html
    assert "has no upstream stop codon" in html


def test_ambiguous_input_carries_the_approximate_mass_note():
    html = render_html(
        [RecordReport(Record("n", "", "ACGTNN"), metrics("ACGTNN"), [])],
        source="s", parameters=PARAMETERS, version="1.0.0",
    )
    assert "Mass is approximate" in html
