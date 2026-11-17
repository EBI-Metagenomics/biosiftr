#!/usr/bin/env python3

# Copyright 2025 EMBL - European Bioinformatics Institute
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

"""Build the BioSIFTR functional profile table of one MGnify species.

Mode ``pre`` writes the protein FASTA of the accessory genes (pangenome clusters that are
absent from the species representative). These genes need an eggNOG annotation. Each cluster
is the protein of one member gene, translated from the CDS in the member GFF.

Mode ``post`` writes ``<rep>_clstr.tsv``: the genes of the representative with their
eggNOG annotation, plus one member gene for each accessory cluster with the eggNOG
annotation of the cluster. Each gene is labelled as core (``true``) or not (``false``).
"""

from __future__ import annotations

import argparse
import csv
import gzip
import logging
import sys
from collections import defaultdict
from pathlib import Path
from typing import IO

from Bio.Seq import Seq

CSV_META_COLUMNS = 3  # Gene, Non-unique Gene name, Annotation
EVALUE_THRESHOLD = 1e-10
TRANSLATION_TABLE = 11  # bacterial and archaeal code (Prokka)
OUTPUT_COLUMNS = ["#contig", "gene_id", "start", "end", "strand", "kegg", "cazy", "pfam", "core"]


def open_text(path: Path) -> IO[str]:
    """Open a text file for reading. The file can be compressed (``.gz``).

    :param path: The file to open.
    :return: A text file handle.
    """
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return open(path)


def split_cell(cell: str) -> list[str]:
    """Split a genome cell of the CSV into locus tags.

    Old mmseqs CSV files join paralogs with a TAB, Panaroo and new mmseqs files use ``;``.

    :param cell: The cell value.
    :return: The locus tags. The list is empty if the cluster is absent.
    """
    return [tag for tag in cell.replace("\t", ";").split(";") if tag]


def pfam_parser(pfam_dat: Path) -> dict[str, str]:
    """Map the Pfam names to the Pfam accessions (without version).

    :param pfam_dat: The ``Pfam-A.hmm.dat.gz`` file.
    :return: {Pfam name: Pfam accession}.
    """
    pfam_desc = {}
    entry = {}
    with open_text(pfam_dat) as handle:
        for line in handle:
            line = line.strip()
            if line == "//":
                if "ID" in entry and "AC" in entry:
                    pfam_desc[entry["ID"]] = entry["AC"]
                entry = {}
            elif line.startswith("#=GF ID"):
                entry["ID"] = line.split()[-1]
            elif line.startswith("#=GF AC"):
                entry["AC"] = line.split()[-1].split(".")[0]
    return pfam_desc


def read_pangenome_csv(path: Path, rep: str) -> tuple[list[str], list[tuple[str, list[str], tuple[str, str] | None]]]:
    """Read the gene_presence_absence.csv of a species one row at a time.

    The genome cells are not kept, because a large species has millions of them.

    :param path: The CSV file. It can be compressed (``.gz``).
    :param rep: The species representative accession.
    :return: The genome columns, and one (cluster, representative locus tags, first member) tuple
        for each row. The first member is the (genome, locus tag) of the first genome column that
        has the cluster. It is set only for the accessory clusters (absent from the representative).
    :raises ValueError: If the representative is not a column of the CSV, or the names are not unique.
    """
    csv.field_size_limit(sys.maxsize)
    clusters = []
    names = set()
    with open_text(path) as handle:
        reader = csv.reader(handle)
        header = next(reader)
        if rep not in header:
            raise ValueError(f"{rep} is not a genome column of {path}")
        rep_index = header.index(rep)
        genomes = header[CSV_META_COLUMNS:]
        for row in reader:
            if row[0] in names:
                raise ValueError(
                    f"{path} has duplicate gene names. Correct the pangenome files first (pangenome_fix.py)"
                )
            names.add(row[0])
            rep_tags = split_cell(row[rep_index])
            clusters.append((row[0], rep_tags, None if rep_tags else first_member(genomes, row)))
    return genomes, clusters


def first_member(genomes: list[str], row: list[str]) -> tuple[str, str] | None:
    """Get the first member gene of a cluster (column order).

    :param genomes: The genome columns of the CSV.
    :param row: The CSV row of the cluster.
    :return: The (genome, locus tag) of the first genome that has the cluster, or ``None``.
    """
    for genome, cell in zip(genomes, row[CSV_META_COLUMNS:]):
        tags = split_cell(cell)
        if tags:
            return genome, tags[0]
    return None


def read_core_clusters(prevalence: Path) -> set[str]:
    """Read the clusters that CGT labels as core.

    :param prevalence: The ``gene_prevalence_corrected.txt`` file (``gene``, ``count``, ``label``).
    :return: The names of the core clusters.
    """
    core = set()
    with open_text(prevalence) as handle:
        for line in handle:
            if line.startswith("#") or line.startswith("gene\t"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3 and parts[2] == "core":
                core.add(parts[0])
    return core


def gff_parser(gff: Path, wanted: set[str] | None = None) -> dict[str, list[str]]:
    """Read the CDS coordinates of a GFF file.

    :param gff: The GFF file. It can be compressed (``.gz``).
    :param wanted: Keep only these locus tags. Keep all CDS if it is ``None``.
    :return: {locus tag: [contig, start, end, strand]}. The strand is ``1`` or ``-1``.
    """
    genes = {}
    with open_text(gff) as handle:
        for line in handle:
            if line.startswith("##FASTA"):
                break
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9 or parts[2] != "CDS":
                continue
            locus_tag = parts[8].split(";")[0].replace("ID=", "")
            if wanted is None or locus_tag in wanted:
                strand = "1" if parts[6] == "+" else "-1"
                genes[locus_tag] = [parts[0], parts[3], parts[4], strand]
    return genes


def eggnog_parser(eggnog: Path, pfam_desc: dict[str, str]) -> dict[str, tuple[str, str, str]]:
    """Read the KEGG, CAZy and Pfam annotations of an eggNOG-mapper output.

    Only the hits with an e-value of 1e-10 or less are kept.

    :param eggnog: The eggNOG-mapper annotations (``*_eggNOG.tsv`` or ``*.emapper.annotations``).
    :param pfam_desc: Output of :func:`pfam_parser`.
    :return: {query: (kegg, cazy, pfam)}.
    """
    annotations = {}
    with open_text(eggnog) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 21:
                continue
            try:
                evalue = float(parts[2])
            except ValueError:  # header line without "#"
                continue
            if evalue > EVALUE_THRESHOLD:
                continue
            kegg = parts[11].replace("ko:", "")
            pfam_ids = sorted({pfam_desc[name] for name in parts[20].split(",") if name in pfam_desc})
            annotations[parts[0]] = (kegg, parts[18], ",".join(pfam_ids) if pfam_ids else "-")
    return annotations


def gff_contigs(gff: Path, wanted: set[str]) -> dict[str, str]:
    """Read the contig sequences of the ``##FASTA`` section of a GFF file.

    :param gff: The GFF file. It can be compressed (``.gz``).
    :param wanted: Keep only these contigs.
    :return: {contig: sequence}.
    """
    sequences = {}
    name, chunks = None, []
    with open_text(gff) as handle:
        for line in handle:
            if line.startswith("##FASTA"):
                break
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    sequences[name] = "".join(chunks)
                name, chunks = line[1:].split()[0], []
                if name not in wanted:
                    name = None  # skip the sequence lines of this record
            elif line and name is not None:
                chunks.append(line)
    if name is not None:
        sequences[name] = "".join(chunks)
    return sequences


def cds_protein(contig: str, start: str, end: str, strand: str) -> str:
    """Translate one CDS.

    :param contig: The contig sequence.
    :param start: The CDS start (1-based, GFF).
    :param end: The CDS end (1-based, inclusive, GFF).
    :param strand: ``1`` or ``-1``.
    :return: The protein sequence, without the final stop codon.
    """
    cds = Seq(contig[int(start) - 1 : int(end)].upper())
    if strand == "-1":
        cds = cds.reverse_complement()
    cds = cds[: len(cds) - len(cds) % 3]  # partial genes at the contig edges
    return str(cds.translate(table=TRANSLATION_TABLE)).rstrip("*")


def write_accessory_fasta(csv_file: Path, gff_dir: Path, rep: str, output: Path) -> int:
    """Write the protein FASTA of the accessory genes (mode ``pre``).

    Each accessory cluster gets the protein of its first member gene (column order). Mode
    ``post`` uses the same member gene.

    :param csv_file: The gene_presence_absence.csv of the species.
    :param gff_dir: The folder with the member GFF files. The GFF files must have the contigs.
    :param rep: The species representative accession.
    :param output: The output FASTA.
    :return: The number of accessory genes in the FASTA.
    """
    _, clusters = read_pangenome_csv(csv_file, rep)
    by_genome = defaultdict(dict)  # genome -> {locus tag: cluster}
    for cluster, _, member in clusters:
        if member is not None:
            genome, locus_tag = member
            by_genome[genome][locus_tag] = cluster

    count = 0
    with open(output, "w") as handle:
        for genome in sorted(by_genome):
            gff = gff_dir / f"{genome}.gff.gz"
            members = by_genome[genome]
            coordinates = gff_parser(gff, set(members))
            contigs = gff_contigs(gff, {coords[0] for coords in coordinates.values()})
            for locus_tag, cluster in members.items():
                coords = coordinates.get(locus_tag)
                if coords is None or coords[0] not in contigs:
                    logging.warning("%s: no CDS sequence for %s (cluster %s)", rep, locus_tag, cluster)
                    continue
                handle.write(f">{cluster}\n{cds_protein(contigs[coords[0]], *coords[1:])}\n")
                count += 1
    return count


def build_profile(
    rep: str,
    rep_gff: Path,
    rep_eggnog: Path,
    pfam_desc: dict[str, str],
    csv_file: Path | None = None,
    prevalence: Path | None = None,
    gff_dir: Path | None = None,
    emapper: Path | None = None,
) -> tuple[int, list[list[str]]]:
    """Build the functional profile rows of one species (mode ``post``).

    A singleton species needs only the representative files. All its genes are core.

    :param rep: The species representative accession.
    :param rep_gff: The GFF of the representative.
    :param rep_eggnog: The eggNOG annotation of the representative.
    :param pfam_desc: Output of :func:`pfam_parser`.
    :param csv_file: The gene_presence_absence.csv (pangenome species only).
    :param prevalence: The CGT ``gene_prevalence_corrected.txt`` (pangenome species only).
    :param gff_dir: The folder with the member GFF files (pangenome species only).
    :param emapper: The eggNOG-mapper annotations of the accessory genes. It can be missing
        when the species has no accessory genes.
    :return: The cluster size and the output rows.
    """
    rep_genes = gff_parser(rep_gff)
    rep_annotations = eggnog_parser(rep_eggnog, pfam_desc)

    def annotated(locus_tag: str, coordinates: list[str], annotation: tuple[str, str, str], core: bool) -> list[str]:
        contig, start, end, strand = coordinates
        return [contig, locus_tag, start, end, strand, *annotation, "true" if core else "false"]

    no_hit = ("-", "-", "-")
    if csv_file is None:
        rows = [annotated(tag, coords, rep_annotations.get(tag, no_hit), True) for tag, coords in rep_genes.items()]
        return 1, rows

    genomes, clusters = read_pangenome_csv(csv_file, rep)
    core_clusters = read_core_clusters(prevalence)
    core_rep_genes = {tag for cluster, rep_tags, _ in clusters if cluster in core_clusters for tag in rep_tags}
    rows = [
        annotated(tag, coords, rep_annotations.get(tag, no_hit), tag in core_rep_genes)
        for tag, coords in rep_genes.items()
    ]

    # One member gene per accessory cluster: the first genome (column order) that has it
    selected = {}  # cluster -> member locus tag
    by_genome = defaultdict(set)
    for cluster, _, member in clusters:
        if member is not None:
            genome, locus_tag = member
            selected[cluster] = locus_tag
            by_genome[genome].add(locus_tag)

    coordinates = {}
    for genome in sorted(by_genome):
        coordinates.update(gff_parser(gff_dir / f"{genome}.gff.gz", by_genome[genome]))

    acc_annotations = eggnog_parser(emapper, pfam_desc) if emapper else {}
    for cluster, locus_tag in selected.items():
        if locus_tag not in coordinates:
            logging.warning("%s: no CDS coordinates for %s (cluster %s)", rep, locus_tag, cluster)
            continue
        annotation = acc_annotations.get(cluster, no_hit)
        rows.append(annotated(locus_tag, coordinates[locus_tag], annotation, cluster in core_clusters))
    return len(genomes), rows


def write_profile(output: Path, cluster_size: int, rows: list[list[str]]) -> None:
    """Write the ``<rep>_clstr.tsv`` table.

    :param output: The output file.
    :param cluster_size: The number of genomes of the species.
    :param rows: Output of :func:`build_profile`.
    """
    with open(output, "w") as handle:
        handle.write(f"#cluster size = {cluster_size}\n")
        handle.write("\t".join(OUTPUT_COLUMNS) + "\n")
        for row in rows:
            handle.write("\t".join(row) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--mode", choices=["pre", "post"], required=True, help="pre: accessory FASTA; post: profile")
    parser.add_argument("--rep", required=True, help="The species representative accession")
    parser.add_argument("--csv", type=Path, help="gene_presence_absence.csv(.gz); not used for singletons")
    parser.add_argument("--prevalence", type=Path, help="CGT gene_prevalence_corrected.txt(.gz) (mode post)")
    parser.add_argument("--rep_gff", type=Path, help="GFF of the representative (mode post)")
    parser.add_argument("--rep_eggnog", type=Path, help="eggNOG annotation of the representative (mode post)")
    parser.add_argument("--gff_dir", type=Path, help="Folder with the member GFF files")
    parser.add_argument("--emapper", type=Path, help="eggNOG-mapper annotations of the accessory genes (mode post)")
    parser.add_argument("--pfam_dat", type=Path, help="Pfam-A.hmm.dat.gz (mode post)")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    if args.mode == "pre":
        if not (args.csv and args.gff_dir):
            parser.error("mode pre needs --csv and --gff_dir")
        count = write_accessory_fasta(args.csv, args.gff_dir, args.rep, Path(f"{args.rep}_accessory.fasta"))
        logging.info("%s: %d accessory genes", args.rep, count)
        return

    if not (args.rep_gff and args.rep_eggnog and args.pfam_dat):
        parser.error("mode post needs --rep_gff, --rep_eggnog and --pfam_dat")
    if args.csv and not (args.prevalence and args.gff_dir):
        parser.error("mode post with --csv needs --prevalence and --gff_dir")
    cluster_size, rows = build_profile(
        args.rep,
        args.rep_gff,
        args.rep_eggnog,
        pfam_parser(args.pfam_dat),
        args.csv,
        args.prevalence,
        args.gff_dir,
        args.emapper,
    )
    write_profile(Path(f"{args.rep}_clstr.tsv"), cluster_size, rows)
    logging.info("%s: %d genes in the profile", args.rep, len(rows))


if __name__ == "__main__":
    main()
