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

"""Tests for bin/pangenome_fix.py."""

from __future__ import annotations

import gzip
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "bin"))

from pangenome_fix import (  # noqa: E402
    CSV_NAME,
    FASTA_NAME,
    RTAB_NAME,
    PangenomeFixError,
    extract_cds,
    find_duplicates,
    fix_pangenome,
    new_names_for_duplicates,
    read_csv,
    read_fasta,
    split_cell,
)

# Contig of genome G1. The CDS of each gene are slices of it.
G1_CONTIG = "ATGAAATTTTAA" + "ATGCCCGGGTAA" + "ATGTTTCCCTAA" + "TTAGGGCCCCAT"
G2_CONTIG = "ATGAAATTTTAA" + "ATGTTTCCCTAA"

GFF_G1 = (
    "##gff-version 3\n"
    "c1\tProkka\tCDS\t1\t12\t.\t+\t0\tID=G1_00001;gene=gyrA\n"
    "c1\tProkka\tCDS\t13\t24\t.\t+\t0\tID=G1_00002;gene=gyrA_2\n"
    "c1\tProkka\tCDS\t25\t36\t.\t+\t0\tID=G1_00003;gene=gyrA\n"
    "c1\tProkka\tCDS\t37\t48\t.\t-\t0\tID=G1_00004;gene=recA\n"
    f"##FASTA\n>c1\n{G1_CONTIG}\n"
)
GFF_G2 = (
    "##gff-version 3\n"
    "c2\tProkka\tCDS\t1\t12\t.\t+\t0\tID=G2_00001;gene=gyrA\n"
    "c2\tProkka\tCDS\t13\t24\t.\t+\t0\tID=G2_00002;gene=gyrA\n"
    f"##FASTA\n>c2\n{G2_CONTIG}\n"
)

# Cluster gyrA (G1_00001, G2_00001) owns "gyrA".
# Cluster with dominant name gyrA (G1_00003, G2_00002) got the suffix "gyrA_2" (bug).
# Cluster with Prokka name gyrA_2 (G1_00002) owns "gyrA_2".
CSV_BAD = (
    "Gene,Non-unique Gene name,Annotation,G1,G2\n"
    "gyrA,gyrA,DNA gyrase,G1_00001,G2_00001\n"
    "gyrA_2,gyrA,DNA gyrase,G1_00003,G2_00002\n"
    "gyrA_2,,DNA gyrase,G1_00002,\n"
    "recA,,Recombinase,G1_00004\tG1_00001,\n"
)
# Records in a different order than the CSV, as in real files.
FASTA_BAD = (
    ">gyrA_2\nATGCCCGGGTAA\n"  # sequence of G1_00002 -> Prokka gyrA_2 row
    ">gyrA\nATGAAATTTTAA\n"
    ">recA\nATGGGGCCCTAA\n"  # reverse complement of the last 12 bases of c1
    ">gyrA_2\nATGTTTCCCTAA\n"  # sequence of G1_00003 -> suffixed row
)
RTAB_BAD = "Gene\tG1\tG2\ngyrA\t1\t1\ngyrA_2\t1\t0\nrecA\t1\t0\n"


def write_gz(path: Path, text: str) -> None:
    with gzip.open(path, "wt") as handle:
        handle.write(text)


@pytest.fixture
def species(tmp_path: Path) -> tuple[Path, Path]:
    pangenome_dir = tmp_path / "pan-genome"
    gff_dir = tmp_path / "genomes1"
    pangenome_dir.mkdir()
    gff_dir.mkdir()
    write_gz(pangenome_dir / CSV_NAME, CSV_BAD)
    write_gz(pangenome_dir / RTAB_NAME, RTAB_BAD)
    write_gz(pangenome_dir / FASTA_NAME, FASTA_BAD)
    write_gz(gff_dir / "G1.gff.gz", GFF_G1)
    write_gz(gff_dir / "G2.gff.gz", GFF_G2)
    return pangenome_dir, gff_dir


def test_split_cell_accepts_tab_and_semicolon() -> None:
    assert split_cell("A_1\tA_2") == ["A_1", "A_2"]
    assert split_cell("A_1;A_2") == ["A_1", "A_2"]
    assert split_cell("") == []


def test_find_duplicates() -> None:
    assert find_duplicates(["a", "b", "a", "c"]) == {"a": [0, 2]}


def test_new_name_skips_names_in_use() -> None:
    rows = [
        ["gyrA", "gyrA", ""],
        ["gyrA_2", "gyrA", ""],
        ["gyrA_2", "", ""],
        ["gyrA_3", "gyrA", ""],
    ]
    duplicates = find_duplicates([row[0] for row in rows])
    assert new_names_for_duplicates(rows, duplicates) == {1: "gyrA_4"}


def test_new_name_fails_without_owner_row() -> None:
    rows = [["gyrA_2", "gyrA", ""], ["gyrA_2", "gyrA", ""]]
    with pytest.raises(PangenomeFixError):
        new_names_for_duplicates(rows, find_duplicates([row[0] for row in rows]))


def test_extract_cds_minus_strand_and_length_filter(tmp_path: Path) -> None:
    gff = tmp_path / "G1.gff.gz"
    write_gz(gff, GFF_G1)
    sequences = extract_cds(gff, {"G1_00004": {12}, "G1_00001": {99}})
    assert sequences == {"G1_00004": "ATGGGGCCCTAA"}


def test_extract_cds_fails_without_fasta_section(tmp_path: Path) -> None:
    gff = tmp_path / "G1.gff"
    gff.write_text(GFF_G1.split("##FASTA")[0])
    with pytest.raises(PangenomeFixError):
        extract_cds(gff, {"G1_00001": {12}})


def test_consistent_species_writes_nothing(tmp_path: Path) -> None:
    pangenome_dir = tmp_path / "pan-genome"
    pangenome_dir.mkdir()
    write_gz(pangenome_dir / CSV_NAME, "Gene,Non-unique Gene name,Annotation,G1\ngyrA,,x,G1_00001\n")
    write_gz(pangenome_dir / RTAB_NAME, "Gene\tG1\ngyrA\t1\n")
    write_gz(pangenome_dir / FASTA_NAME, ">gyrA\nATG\n")
    outdir = tmp_path / "fixed"
    report = fix_pangenome(pangenome_dir, tmp_path, "MGYG1", outdir)
    assert report["status"] == "ok"
    assert not outdir.exists()


def test_fix_pangenome_end_to_end(species: tuple[Path, Path], tmp_path: Path) -> None:
    pangenome_dir, gff_dir = species
    outdir = tmp_path / "fixed"
    report = fix_pangenome(pangenome_dir, gff_dir, "MGYG1", outdir)
    assert report["status"] == "fixed"
    assert report["renamed_rows"] == "1"

    _, rows = read_csv(outdir / CSV_NAME)
    assert [row[0] for row in rows] == ["gyrA", "gyrA_3", "gyrA_2", "recA"]
    assert rows[3][3] == "G1_00004;G1_00001"

    records = dict(read_fasta(outdir / FASTA_NAME))
    assert records == {
        "gyrA_2": "ATGCCCGGGTAA",
        "gyrA": "ATGAAATTTTAA",
        "recA": "ATGGGGCCCTAA",
        "gyrA_3": "ATGTTTCCCTAA",
    }

    with gzip.open(outdir / RTAB_NAME, "rt") as handle:
        rtab = handle.read().splitlines()
    assert rtab == ["Gene\tG1\tG2", "gyrA\t1\t1", "gyrA_2\t1\t0", "gyrA_3\t1\t1", "recA\t1\t0"]

    renamed = (outdir / "renamed_genes.tsv").read_text().splitlines()
    assert renamed[1] == "gyrA_2\tgyrA_3\tgyrA"


def test_fix_pangenome_output_is_reproducible(species: tuple[Path, Path], tmp_path: Path) -> None:
    pangenome_dir, gff_dir = species
    fix_pangenome(pangenome_dir, gff_dir, "MGYG1", tmp_path / "a")
    fix_pangenome(pangenome_dir, gff_dir, "MGYG1", tmp_path / "b")
    for name in (CSV_NAME, RTAB_NAME, FASTA_NAME):
        assert (tmp_path / "a" / name).read_bytes() == (tmp_path / "b" / name).read_bytes()


def test_fix_pangenome_fails_when_no_member_matches(species: tuple[Path, Path], tmp_path: Path) -> None:
    pangenome_dir, gff_dir = species
    write_gz(pangenome_dir / FASTA_NAME, FASTA_BAD.replace("ATGTTTCCCTAA", "ATGTTTTTTTAA"))
    with pytest.raises(PangenomeFixError):
        fix_pangenome(pangenome_dir, gff_dir, "MGYG1", tmp_path / "fixed")


def test_fix_pangenome_versioned_genomes(species: tuple[Path, Path], tmp_path: Path) -> None:
    # Genome G1.1 has the locus tags G1_<number>: the GFF name comes from the CSV column, not from the tag
    pangenome_dir, gff_dir = species
    write_gz(pangenome_dir / CSV_NAME, CSV_BAD.replace(",G1,G2\n", ",G1.1,G2.1\n"))
    (gff_dir / "G1.gff.gz").rename(gff_dir / "G1.1.gff.gz")
    (gff_dir / "G2.gff.gz").rename(gff_dir / "G2.1.gff.gz")
    report = fix_pangenome(pangenome_dir, gff_dir, "G1.1", tmp_path / "fixed")
    assert report["status"] == "fixed"
    assert dict(read_fasta(tmp_path / "fixed" / FASTA_NAME))["gyrA_3"] == "ATGTTTCCCTAA"


def test_fix_pangenome_fails_when_member_gff_missing(species: tuple[Path, Path], tmp_path: Path) -> None:
    pangenome_dir, gff_dir = species
    (gff_dir / "G2.gff.gz").unlink()
    with pytest.raises(PangenomeFixError):
        fix_pangenome(pangenome_dir, gff_dir, "MGYG1", tmp_path / "fixed")
