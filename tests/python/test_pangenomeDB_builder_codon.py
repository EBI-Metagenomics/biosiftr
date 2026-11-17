# Copyright 2026 EMBL - European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Tests for bin/pangenomeDB_builder_codon.py."""

from __future__ import annotations

import gzip
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "bin"))

from pangenomeDB_builder_codon import (  # noqa: E402
    build_profile,
    pfam_parser,
    write_accessory_fasta,
    write_profile,
)

REP = "R1"
REP_GFF = (
    "##gff-version 3\n"
    "c1\tProkka\tCDS\t1\t300\t.\t+\t0\tID=R1_00001;gene=a\n"
    "c1\tProkka\tCDS\t301\t600\t.\t-\t0\tID=R1_00002;gene=b\n"
    "##FASTA\n>c1\nACGT\n"
)
# Member contig: CDS M2_00002 (100-190, minus) and M2_00004 (300-390, plus). Each CDS is
# 91 bp: start codon, 28 codons, stop codon and one extra base (partial codon).
CDS_C = "ATG" + "GGG" * 28 + "TAA" + "A"
CDS_D = "ATG" + "TTT" * 28 + "TAA" + "A"
COMPLEMENT = str.maketrans("ACGT", "TGCA")
MEMBER_CONTIG = "N" * 99 + CDS_C.translate(COMPLEMENT)[::-1] + "N" * 109 + CDS_D
MEMBER_GFF = (
    "##gff-version 3\n"
    "m1\tProkka\tCDS\t1\t90\t.\t+\t0\tID=M2_00001\n"
    "m1\tProkka\tCDS\t100\t190\t.\t-\t0\tID=M2_00002\n"
    "m1\tProkka\tCDS\t200\t290\t.\t+\t0\tID=M2_00003\n"
    "m1\tProkka\tCDS\t300\t390\t.\t+\t0\tID=M2_00004\n"
    f"##FASTA\n>m0\nACGT\n>m1\n{MEMBER_CONTIG[:200]}\n{MEMBER_CONTIG[200:]}\n"
)
# Cluster c has two paralogs in M2, joined with a TAB as in old mmseqs CSV files.
CSV = (
    "Gene,Non-unique Gene name,Annotation,R1,M2\n"
    "a,,x,R1_00001,M2_00001\n"
    "b,,x,R1_00002,\n"
    "c,,x,,M2_00002\tM2_00003\n"
    "d,,x,,M2_00004\n"
)
PREVALENCE = "# Core threshold: >= 2 observations\ngene\tcount\tlabel\na\t2\tcore\nb\t1\trare\nc\t1\trare\nd\t1\tcore\n"
PFAM_DAT = "# STOCKHOLM 1.0\n#=GF ID   Glyco_hydro_1\n#=GF AC   PF00232.21\n//\n"


def eggnog_line(query: str, evalue: str, kegg: str, cazy: str, pfams: str) -> str:
    """Make one eggNOG-mapper annotation line (21 columns)."""
    columns = [query, "seed", evalue] + ["-"] * 18
    columns[11], columns[18], columns[20] = kegg, cazy, pfams
    return "\t".join(columns) + "\n"


EGGNOG_HEADER = "#query\tseed_ortholog\tevalue\n"
REP_EGGNOG = (
    EGGNOG_HEADER
    + eggnog_line("R1_00001", "1e-50", "ko:K00001", "GH1", "Glyco_hydro_1,Unknown")
    + eggnog_line("R1_00002", "1e-3", "ko:K00009", "-", "-")
)
EMAPPER = EGGNOG_HEADER + eggnog_line("c", "1e-20", "ko:K00002,ko:K00003", "-", "-")


def write_gz(path: Path, text: str) -> None:
    with gzip.open(path, "wt") as handle:
        handle.write(text)


@pytest.fixture
def species(tmp_path: Path) -> dict[str, Path]:
    gff_dir = tmp_path / "genomes1"
    gff_dir.mkdir()
    files = {
        "csv": tmp_path / "gene_presence_absence.csv.gz",
        "prevalence": tmp_path / "gene_prevalence_corrected.txt.gz",
        "rep_gff": tmp_path / "R1.gff",
        "rep_eggnog": tmp_path / "R1_eggNOG.tsv",
        "emapper": tmp_path / "R1_out.emapper.annotations",
        "pfam_dat": tmp_path / "Pfam-A.hmm.dat.gz",
        "gff_dir": gff_dir,
    }
    write_gz(files["csv"], CSV)
    write_gz(files["prevalence"], PREVALENCE)
    files["rep_gff"].write_text(REP_GFF)
    files["rep_eggnog"].write_text(REP_EGGNOG)
    files["emapper"].write_text(EMAPPER)
    write_gz(files["pfam_dat"], PFAM_DAT)
    write_gz(gff_dir / "M2.gff.gz", MEMBER_GFF)
    return files


def test_pre_writes_accessory_proteins(species: dict[str, Path], tmp_path: Path) -> None:
    output = tmp_path / "R1_accessory.fasta"
    assert write_accessory_fasta(species["csv"], species["gff_dir"], REP, output) == 2
    assert output.read_text() == f">c\nM{'G' * 28}\n>d\nM{'F' * 28}\n"


def test_pre_skips_genes_without_cds(species: dict[str, Path], tmp_path: Path) -> None:
    write_gz(species["gff_dir"] / "M2.gff.gz", MEMBER_GFF.replace("ID=M2_00004", "ID=other"))
    output = tmp_path / "R1_accessory.fasta"
    assert write_accessory_fasta(species["csv"], species["gff_dir"], REP, output) == 1
    assert output.read_text().startswith(">c\n")


def test_pre_fails_on_duplicate_names(species: dict[str, Path], tmp_path: Path) -> None:
    write_gz(species["csv"], CSV.replace("d,,x", "c,,x"))
    with pytest.raises(ValueError):
        write_accessory_fasta(species["csv"], species["gff_dir"], REP, tmp_path / "out.fasta")


def test_post_pangenome_species(species: dict[str, Path]) -> None:
    cluster_size, rows = build_profile(
        REP,
        species["rep_gff"],
        species["rep_eggnog"],
        pfam_parser(species["pfam_dat"]),
        species["csv"],
        species["prevalence"],
        species["gff_dir"],
        species["emapper"],
    )
    assert cluster_size == 2
    assert rows == [
        ["c1", "R1_00001", "1", "300", "1", "K00001", "GH1", "PF00232", "true"],
        ["c1", "R1_00002", "301", "600", "-1", "-", "-", "-", "false"],
        ["m1", "M2_00002", "100", "190", "-1", "K00002,K00003", "-", "-", "false"],
        ["m1", "M2_00004", "300", "390", "1", "-", "-", "-", "true"],
    ]


def test_post_singleton_species(species: dict[str, Path], tmp_path: Path) -> None:
    cluster_size, rows = build_profile(REP, species["rep_gff"], species["rep_eggnog"], pfam_parser(species["pfam_dat"]))
    assert cluster_size == 1
    assert [row[-1] for row in rows] == ["true", "true"]

    output = tmp_path / "R1_clstr.tsv"
    write_profile(output, cluster_size, rows)
    lines = output.read_text().splitlines()
    assert lines[0] == "#cluster size = 1"
    assert lines[1].split("\t")[0] == "#contig"
    assert len(lines) == 4


def test_post_without_emapper_output(species: dict[str, Path]) -> None:
    _, rows = build_profile(
        REP,
        species["rep_gff"],
        species["rep_eggnog"],
        pfam_parser(species["pfam_dat"]),
        species["csv"],
        species["prevalence"],
        species["gff_dir"],
    )
    assert [row[1] for row in rows] == ["R1_00001", "R1_00002", "M2_00002", "M2_00004"]
    assert rows[2][5:8] == ["-", "-", "-"]
