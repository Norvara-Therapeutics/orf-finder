import pytest

from orf_finder.genetic_code import (
    CODONS,
    available_tables,
    get_genetic_code,
    resolve_table_name,
)
from orf_finder.translate import CODON_TABLE, translate


def test_every_bundled_table_is_complete_and_well_formed():
    tables = available_tables()
    assert len(tables) == 27
    for code in tables:
        assert len(code.codons) == 64
        assert set(code.codons) == set(CODONS)
        assert code.stop_codons, f"table {code.id} has no stop codon"
        assert code.start_codons, f"table {code.id} has no start codon"
        assert not code.start_codons & code.stop_codons
        assert code.names and code.name == code.names[0]


def test_table_ids_are_the_numbers_ncbi_assigns():
    ids = [code.id for code in available_tables()]
    assert ids == sorted(ids)
    assert ids[0] == 1 and ids[-1] == 33
    # NCBI leaves several numbers unused.
    assert 7 not in ids and 8 not in ids and 17 not in ids


def test_table_one_is_the_standard_code():
    standard = get_genetic_code(1)
    assert dict(standard.codons) == CODON_TABLE
    assert standard.name == "Standard"
    assert standard.stop_codons == {"TAA", "TAG", "TGA"}
    assert standard.start_codons == {"ATG", "CTG", "TTG"}


def test_every_initiation_codon_initiates_as_methionine():
    # The amino-acid marks in NCBI's sncbieaa strings are all 'M'.
    for code in available_tables():
        for codon in code.start_codons:
            assert codon in code.codons


@pytest.mark.parametrize(
    ("table_id", "name", "stops"),
    [
        (27, "Karyorelict Nuclear", {"TGA"}),
        (28, "Condylostoma Nuclear", {"TAA", "TAG", "TGA"}),
        (31, "Blastocrithidia Nuclear", {"TAA", "TAG"}),
    ],
)
def test_readthrough_tables_take_their_stops_from_the_start_string(table_id, name, stops):
    """Tables 27, 28 and 31 translate every stop triplet as a sense codon.

    Termination in those organisms is context-dependent, so NCBI records the
    potential terminator in sncbieaa rather than ncbieaa. Reading only ncbieaa
    would leave the table with no stop codon and no ORF could ever terminate.
    """
    code = get_genetic_code(table_id)
    assert code.name == name
    assert code.stop_codons == stops
    for codon in stops:
        assert code.codons[codon] != "*", "the codon still translates as a sense codon"


def test_alternative_tables_reassign_codons():
    mito = get_genetic_code(2)
    assert mito.codons["TGA"] == "W"
    assert mito.stop_codons == {"TAA", "TAG", "AGA", "AGG"}
    ciliate = get_genetic_code(6)
    assert ciliate.codons["TAA"] == "Q" and ciliate.codons["TAG"] == "Q"
    assert ciliate.stop_codons == {"TGA"}


def test_bacterial_table_translates_like_the_standard_code_but_starts_differ():
    bacterial = get_genetic_code(11)
    assert dict(bacterial.codons) == CODON_TABLE
    assert {"GTG", "TTG", "ATT", "ATC", "ATA"} <= bacterial.start_codons
    assert "GTG" not in get_genetic_code(1).start_codons


@pytest.mark.parametrize(
    ("value", "table_id"),
    [
        ("1", 1), ("11", 11), ("standard", 1), ("Standard", 1), ("bacterial", 11),
        ("plant-plastid", 11), ("vertebrate-mito", 2),
        ("Vertebrate Mitochondrial", 2), ("vertebrate_mitochondrial", 2),
        ("ciliate-nuclear", 6),
    ],
)
def test_tables_resolve_by_number_name_and_alias(value, table_id):
    assert resolve_table_name(value) == table_id
    assert get_genetic_code(value).id == table_id


@pytest.mark.parametrize("value", ["0", "7", "99", "martian", ""])
def test_unknown_tables_raise(value):
    with pytest.raises(KeyError):
        get_genetic_code(value)


def test_translation_honours_the_chosen_table():
    assert translate("TGA") == "*"
    assert translate("TGA", get_genetic_code(2)) == "W"
    assert translate("ATGTGATAA", get_genetic_code(2)) == "MW*"
