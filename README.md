# orf-finder

Nucleic-acid metrics and six-frame open reading frame scanning from the command
line. Pure Python with no third-party dependencies, and a self-contained HTML
report that opens from a local disk without a network connection.

## Install

```bash
pip install .
```

## Quickstart

```bash
orf-finder scan --input examples/phix174.fasta --min-length 300 --report phix174.html
```

```
> NC_001422.1  Escherichia phage phiX174, complete genome
  5,386 nt DNA  A=1,291 T=1,684 C=1,157 G=1,254 N=0  GC=44.76%  mass=1,663,975.55 g/mol
  8 ORFs
    ORF      FRAME       START         END       NT      AA  START STOP  
    ORF1        +2       1,001       2,284    1,284     427  ATG   TGA   
    ORF2        +3       2,931       3,917      987     328  ATG   TAA   
    ORF3        +1       3,076       3,684      609     202  ATG   TGA   
    ORF4        -1       1,094       1,651      558     185  ATG   TGA   
    ORF5        +1       2,395       2,922      528     175  ATG   TGA   
    ORF6        +3         390         848      459     152  ATG   TAA   
    ORF7        -1       1,772       2,167      396     131  ATG   TAA   
    ORF8        +1       3,946       4,263      318     105  ATG   TAA   
```

The report adds the full composition breakdown, a per-frame ORF census, the
coordinate table and every translated peptide.

## Options

```
-i, --input PATH         FASTA file to scan, or '-' for standard input
-m, --min-length N       smallest ORF to report, in nucleotides
                         including the stop codon (default: 300)
-g, --table TABLE        genetic code by number, name or alias
                         (default: 1, the standard code)
-b, --begin N            first base of the fragment to scan (default: 1)
-e, --end N              last base of the fragment; 0 means to the end
    --strand {both,plus,minus}
    --start-codons       atg | alt | any | an explicit comma-separated list
-n, --ignore-nested      drop ORFs completely placed within a longer ORF
                         on the same strand
    --include-truncated  also report an ORF that reaches the end of the
                         sequence without a stop codon
    --no-peptides        omit peptide sequences from the HTML report
```

`orf-finder tables` lists the 27 bundled genetic codes.

## Outputs

| Flag | Output |
| --- | --- |
| `--format table` (default) | metrics and a readable ORF table on stdout |
| `--format tsv` \| `json` \| `none` | the same on stdout, or nothing |
| `--output PATH` | the ORF table as TSV |
| `--protein-fasta PATH` | translated peptides as FASTA |
| `--cds-fasta PATH` | ORF coding sequences as nucleotide FASTA |
| `--feature-table PATH` | a five-column feature table |
| `--report PATH` | the self-contained HTML report |

Any combination may be requested in a single pass; the sequence is scanned once.
Feature-table output requires unique identifiers in the input FASTA, so its
sequence names match the source. For a nonstandard genetic code, each CDS carries
the selected `transl_table` qualifier.

## Reading the output

- Coordinates are 1-based and inclusive on the forward strand, with
  `start <= end` on both strands; a minus-strand ORF is read from `end` towards
  `start`.
- Frames `+1`, `+2`, `+3` are offsets into the sequence as given; `-1`, `-2`,
  `-3` are offsets into its reverse complement.
- `length_nt` includes the stop codon and `--min-length` filters on it;
  `length_aa` is the peptide, which excludes it.
- With `--start-codons atg` or `alt`, one ORF is reported per in-frame stop
  codon, beginning at the 5'-most start codon after the preceding stop, and an
  initiation codon is translated as methionine. With `any`, every region between
  two stops is an ORF and the first region of a frame is marked
  `five_prime_partial`, having no upstream stop. These regions have no required
  initiation codon; `has_initiation_codon` is false and feature-table/CDS FASTA
  outputs mark their 5' ends as partial.
- An explicitly listed start codon that the selected genetic code does not
  recognise as an initiator is still scanned, but retains its ordinary amino
  acid and is marked 5' partial.
- `truncated` marks an ORF that reaches the end of the sequence without a stop
  codon.
- GC content is computed over unambiguous bases; mass is approximate when
  ambiguity codes are present. RNA input is detected when U is present and T is
  absent, and is translated through the same table.

## Genetic codes

All 27 translation tables are bundled, transcribed verbatim from the `ncbieaa`
and `sncbieaa` strings of
[gc.prt](https://ftp.ncbi.nlm.nih.gov/entrez/misc/data/gc.prt), whose sha256 is
recorded in `genetic_code.py`. Each table carries its own stop codons and its own
initiation codons, so `--start-codons alt` means something different under each
one:

```bash
orf-finder scan -i genome.fasta --table 11             # bacterial
orf-finder scan -i mito.fasta  --table vertebrate-mito # table 2
```

## Limitations

- Sequences are treated as linear, so a feature crossing the origin of a
  circular molecule is not found.
- Alignment gaps are accepted but not removed, so a gapped sequence yields
  coordinates in alignment space and `X` residues where the gaps fall. Gaps
  contribute to length but are excluded from base composition and mass.
- The whole sequence is held in memory, which suits genes, plasmids, viruses and
  bacterial chromosomes rather than vertebrate chromosomes.
- An ORF is a reading frame, not a gene. No coding potential, ribosome binding
  site, splice site or homology evidence is considered.

## Development

```bash
pip install -e ".[dev]"
pytest -q
python -m pyflakes src/orf_finder tests
```

## License

Licensed under the
[PolyForm Noncommercial License 1.0.0](https://polyformproject.org/licenses/noncommercial/1.0.0).
See [LICENSE](LICENSE) for the complete terms.
