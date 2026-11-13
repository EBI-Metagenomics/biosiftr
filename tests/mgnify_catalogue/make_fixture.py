#!/usr/bin/env python3

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

"""Make the tiny MGnify catalogue used by tests/build_mgnify_db.nf.test.

The catalogue has the same layout as a catalogue folder on the cluster
(``<out>/test-host/v1.0/ftp``) and three species:

- MGYG000000101.1: a singleton with a versioned accession (folder prefix MGYG0000001).
- MGYG000000201: a pangenome with consistent files. Its accessory clusters are the
  SARS-CoV-2 S gene (two paralogs in one genome) and ORF3a.
- MGYG000000301: a pangenome with the duplicate-name bug of old mmseqs catalogues.
  The suffixed cluster (named gyrA_2, like the Prokka gyrA_2 cluster) is the SARS-CoV-2 N gene.

The SARS-CoV-2 genes are used because the tiny eggNOG database (tests/eggnog_db) has
annotations for them. The other genes are random sequences.

Usage::

    make_fixture.py --genome genome.fasta --gff genome.gff3 --out tests/mgnify_catalogue

The SARS-CoV-2 files come from nf-core/test-datasets (modules branch,
data/genomics/sarscov2/genome/).
"""

from __future__ import annotations

import argparse
import gzip
import io
import random
from pathlib import Path

COMPLEMENT = str.maketrans("ACGT", "TGCA")
METADATA_COLUMNS = [
    "Genome",
    "Genome_type",
    "Length",
    "N_contigs",
    "N50",
    "GC_content",
    "Completeness",
    "Contamination",
    "rRNA_5S",
    "rRNA_16S",
    "rRNA_23S",
    "tRNAs",
    "Genome_accession",
    "Species_rep",
    "Lineage",
    "Sample_accession",
    "Study_accession",
    "Country",
    "Continent",
    "FTP_download",
]
EGGNOG_HEADER = [
    "#query",
    "seed_ortholog",
    "evalue",
    "score",
    "eggNOG_OGs",
    "max_annot_lvl",
    "COG_category",
    "Description",
    "Preferred_name",
    "GOs",
    "EC",
    "KEGG_ko",
    "KEGG_Pathway",
    "KEGG_Module",
    "KEGG_Reaction",
    "KEGG_rclass",
    "BRITE",
    "KEGG_TC",
    "CAZy",
    "BiGG_Reaction",
    "PFAMs",
]


def write_gz(path: Path, text: str) -> None:
    """Write a compressed text file with no time stamp (same content, same bytes)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.GzipFile(path, "wb", mtime=0) as handle:
        handle.write(text.encode())


def write_text(path: Path, text: str) -> None:
    """Write a text file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def sarscov2_cds(genome: Path, gff: Path) -> dict[str, str]:
    """Get the CDS sequences of the SARS-CoV-2 genes (first CDS of each gene)."""
    sequence = "".join(line.strip() for line in genome.read_text().splitlines() if not line.startswith(">"))
    cds = {}
    for line in gff.read_text().splitlines():
        parts = line.split("\t")
        if len(parts) == 9 and parts[2] == "CDS":
            gene = dict(item.split("=", 1) for item in parts[8].split(";") if "=" in item)["gene"]
            cds.setdefault(gene, sequence[int(parts[3]) - 1 : int(parts[4])])
    return cds


def random_cds(length: int) -> str:
    """Make a random CDS with a start and a stop codon."""
    body = "".join(random.choice("ACGT") for _ in range(length - 6))
    return "ATG" + body.replace("TAA", "TAC").replace("TAG", "TAC").replace("TGA", "TAC") + "TAA"


def genome_files(accession: str, genes: list[tuple[str, str, str, str]]) -> tuple[str, str]:
    """Make the contig and the Prokka-style GFF (with ##FASTA) of a genome.

    :param genes: (locus tag, gene name, CDS, strand) for each gene.
    :return: The contig sequence and the GFF text.
    """
    contig_name = f"{accession}_1"
    contig, lines = "", []
    for locus_tag, gene_name, cds, strand in genes:
        piece = cds if strand == "+" else cds.translate(COMPLEMENT)[::-1]
        start = len(contig) + 11
        contig += "N" * 10 + piece
        attributes = f"ID={locus_tag}" + (f";gene={gene_name}" if gene_name else "") + ";product=protein"
        lines.append(f"{contig_name}\tProkka\tCDS\t{start}\t{len(contig)}\t.\t{strand}\t0\t{attributes}")
    contig += "N" * 10
    gff = "##gff-version 3\n" + "\n".join(lines) + f"\n##FASTA\n>{contig_name}\n{contig}\n"
    return contig, gff


def eggnog_table(annotations: list[tuple[str, str, str, str]]) -> str:
    """Make an eggNOG-mapper table: (query, KEGG_ko, CAZy, PFAMs) for each annotated gene."""
    lines = ["\t".join(EGGNOG_HEADER)]
    for query, kegg, cazy, pfams in annotations:
        columns = [query, "seed", "1e-50"] + ["-"] * 18
        columns[11], columns[18], columns[20] = kegg, cazy, pfams
        lines.append("\t".join(columns))
    return "\n".join(lines) + "\n"


class Catalogue:
    """Write the species of a tiny catalogue."""

    def __init__(self, root: Path) -> None:
        self.root = root
        self.metadata = ["\t".join(METADATA_COLUMNS)]

    def add_metadata(self, genome: str, rep: str, completeness: float) -> None:
        values = [genome, "MAG", "1000", "1", "1000", "50", str(completeness), "0.5", "0", "0", "0", "20"]
        values += [genome, rep, "d__Bacteria;p__;c__;o__;f__;g__;s__", "SAMN0", "PRJN0", "UK", "Europe", "-"]
        self.metadata.append("\t".join(values))

    def add_species(
        self,
        rep: str,
        members: list[str],
        clusters: list[tuple[str, str, str, dict[str, list[str]]]],
        rep_annotations: list[tuple[str, str, str, str]],
    ) -> None:
        """Add a species.

        :param clusters: (Gene, Non-unique Gene name, CDS, {genome: [locus tags]}) for each cluster.
            The Prokka gene name of the members is the Gene name without a suffix added by the bug.
        :param rep_annotations: The eggNOG annotations of the representative genes.
        """
        genomes = [rep] + members
        prefix = rep.split(".")[0][:-2]  # MGYG000000101.1 -> MGYG0000001
        species_dir = self.root / "species_catalogue" / prefix / rep
        for genome in genomes:
            self.add_metadata(genome, rep, 95 - 5 * genomes.index(genome))

        genes = {genome: [] for genome in genomes}
        for index, (gene, non_unique, cds, members_of) in enumerate(clusters):
            prokka_name = non_unique or gene
            for genome, locus_tags in members_of.items():
                for locus_tag in locus_tags:
                    strand = "-" if index % 2 else "+"
                    genes[genome].append((locus_tag, "" if gene.startswith("group_") else prokka_name, cds, strand))

        for genome in genomes:
            contig, gff = genome_files(genome, sorted(genes[genome]))
            if genome == rep:
                write_text(species_dir / "genome" / f"{rep}.fna", f">{rep}_1\n{contig}\n")
                write_text(species_dir / "genome" / f"{rep}.gff", gff)
                write_text(species_dir / "genome" / f"{rep}_eggNOG.tsv", eggnog_table(rep_annotations))
            if members:
                write_gz(self.root / "all_genomes" / prefix / rep / "genomes1" / f"{genome}.gff.gz", gff)
        if not members:
            return

        # Pangenome files as the old mmseqs code wrote them: paralogs joined with a TAB,
        # Rtab rows keyed by name (a duplicate name keeps only the last cluster), FASTA in another order
        csv_lines = [",".join(["Gene", "Non-unique Gene name", "Annotation"] + genomes)]
        rtab = {}
        for gene, non_unique, _cds, members_of in clusters:
            cells = ["\t".join(members_of.get(genome, [])) for genome in genomes]
            csv_lines.append(",".join([gene, non_unique, "protein"] + cells))
            rtab[gene] = ["1" if members_of.get(genome) else "0" for genome in genomes]
        fasta = [(gene, cds) for gene, _non_unique, cds, _members in reversed(clusters)]
        prevalence = "# Core threshold: >= 2 observations or >= 100.00% frequency\n"
        prevalence += "# Rare threshold: <= 0 observations or <= 0.00% frequency\ngene\tcount\tlabel\n"
        for gene, values in rtab.items():
            label = "core" if all(value == "1" for value in values) else "middle"
            prevalence += f"{gene}\t{values.count('1')}\t{label}\n"

        pangenome_dir = species_dir / "pan-genome"
        write_gz(pangenome_dir / "gene_presence_absence.csv.gz", "\n".join(csv_lines) + "\n")
        rtab_text = io.StringIO()
        rtab_text.write("\t".join(["Gene"] + genomes) + "\n")
        for gene in sorted(rtab):
            rtab_text.write("\t".join([gene] + rtab[gene]) + "\n")
        write_gz(pangenome_dir / "gene_presence_absence.Rtab.gz", rtab_text.getvalue())
        write_gz(pangenome_dir / "pan-genome.fna.gz", "".join(f">{gene}\n{cds}\n" for gene, cds in fasta))
        write_gz(pangenome_dir / "gene_prevalence_corrected.txt.gz", prevalence)

    def write_metadata(self) -> None:
        write_text(self.root / "genomes-all_metadata.tsv", "\n".join(self.metadata) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--genome", type=Path, required=True, help="SARS-CoV-2 genome.fasta")
    parser.add_argument("--gff", type=Path, required=True, help="SARS-CoV-2 genome.gff3")
    parser.add_argument("--out", type=Path, required=True, help="Output folder (tests/mgnify_catalogue)")
    args = parser.parse_args()

    random.seed(42)
    sars = sarscov2_cds(args.genome, args.gff)
    catalogue = Catalogue(args.out / "test-host" / "v1.0" / "ftp")
    kegg_rep = [("ko:K00134", "GH1", "Glyco_hydro_1")]

    # Singleton
    rep = "MGYG000000101.1"
    catalogue.add_species(
        rep,
        [],
        [("pgk", "", random_cds(300), {rep: [f"{rep}_00001"]}), ("tpiA", "", random_cds(240), {rep: [f"{rep}_00002"]})],
        [(f"{rep}_00001", *kegg_rep[0])],
    )

    # Consistent pangenome
    rep, member = "MGYG000000201", "MGYG000000202"
    catalogue.add_species(
        rep,
        [member],
        [
            ("recA", "", random_cds(270), {rep: [f"{rep}_00001"], member: [f"{member}_00001"]}),
            ("dnaA", "", random_cds(210), {rep: [f"{rep}_00002"]}),
            ("group_1", "", sars["S"], {member: [f"{member}_00002", f"{member}_00003"]}),
            ("group_2", "", sars["ORF3a"], {member: [f"{member}_00004"]}),
        ],
        [(f"{rep}_00001", *kegg_rep[0])],
    )

    # Pangenome with the duplicate-name bug
    rep, member = "MGYG000000301", "MGYG000000302"
    catalogue.add_species(
        rep,
        [member],
        [
            ("gyrA", "gyrA", random_cds(270), {rep: [f"{rep}_00001"], member: [f"{member}_00001"]}),
            ("gyrA_2", "", random_cds(180), {rep: [f"{rep}_00002"]}),
            ("gyrA_2", "gyrA", sars["N"], {member: [f"{member}_00002"]}),
        ],
        [(f"{rep}_00001", *kegg_rep[0])],
    )
    catalogue.write_metadata()


if __name__ == "__main__":
    main()
