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

"""Check and correct the mmseqs pangenome files of one MGnify species.

Older releases of the genomes-catalogue-pipeline can give two clusters the same name
(for example a suffixed ``gyrA_2`` and a Prokka paralog name ``gyrA_2``). The CSV and the
pan-genome FASTA then contain the same name twice, and the Rtab loses one of the clusters.

This script:

1. Checks that the CSV, the Rtab and the FASTA have the same clusters and unique names.
2. If they do not, gives the suffixed cluster of each pair a new free name, pairs each
   duplicate FASTA record with its CSV row (by an exact match to a member CDS), and writes
   a corrected CSV, FASTA and an Rtab rebuilt from the CSV.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import logging
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import IO

CSV_NAME = "gene_presence_absence.csv.gz"
RTAB_NAME = "gene_presence_absence.Rtab.gz"
FASTA_NAME = "pan-genome.fna.gz"
CSV_META_COLUMNS = 3  # Gene, Non-unique Gene name, Annotation

COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


class PangenomeFixError(Exception):
    """The pangenome files cannot be corrected without a guess."""


def open_text(path: Path) -> IO[str]:
    """Open a text file for reading. The file can be compressed (``.gz``).

    :param path: The file to open.
    :return: A text file handle.
    """
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return open(path)


def open_gz_writer(path: Path) -> IO[str]:
    """Open a compressed text file for writing.

    The gzip header has no time stamp, so the same content always gives the same file.

    :param path: The file to write.
    :return: A text file handle.
    """
    return io.TextIOWrapper(gzip.GzipFile(path, "wb", mtime=0), newline="")


def read_csv(path: Path) -> tuple[list[str], list[list[str]]]:
    """Read a Roary-style gene_presence_absence.csv file.

    :param path: The CSV file. It can be compressed (``.gz``).
    :return: The header and the rows.
    """
    csv.field_size_limit(sys.maxsize)
    with open_text(path) as handle:
        reader = csv.reader(handle)
        header = next(reader)
        rows = list(reader)
    return header, rows


def count_rtab_rows(path: Path) -> int:
    """Count the gene rows of an Rtab file (all lines except the header).

    :param path: The Rtab file. It can be compressed (``.gz``).
    :return: The number of gene rows.
    """
    with open_text(path) as handle:
        return sum(1 for _ in handle) - 1


def read_fasta(path: Path) -> list[tuple[str, str]]:
    """Read a FASTA file and keep the record order.

    :param path: The FASTA file. It can be compressed (``.gz``).
    :return: A list of (name, sequence). The name is the first word of the header.
    """
    records = []
    name, chunks = None, []
    with open_text(path) as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(chunks)))
                name, chunks = line[1:].split()[0], []
            elif line:
                chunks.append(line)
    if name is not None:
        records.append((name, "".join(chunks)))
    return records


def split_cell(cell: str) -> list[str]:
    """Split a genome cell of the CSV into locus tags.

    Old mmseqs CSV files join paralogs with a TAB, Panaroo and new mmseqs files use ``;``.

    :param cell: The cell value.
    :return: The locus tags. The list is empty if the cluster is absent.
    """
    return [tag for tag in cell.replace("\t", ";").split(";") if tag]


def genome_of(locus_tag: str) -> str:
    """Get the genome accession of a Prokka locus tag (``<genome>_<number>``).

    :param locus_tag: The locus tag.
    :return: The genome accession.
    """
    return locus_tag.rsplit("_", 1)[0]


def find_duplicates(names: list[str]) -> dict[str, list[int]]:
    """Find the names that occur more than one time.

    :param names: The names, in file order.
    :return: {name: [index, ...]} for each name that occurs more than one time.
    """
    positions = defaultdict(list)
    for index, name in enumerate(names):
        positions[name].append(index)
    return {name: indexes for name, indexes in positions.items() if len(indexes) > 1}


def new_names_for_duplicates(rows: list[list[str]], duplicates: dict[str, list[int]]) -> dict[int, str]:
    """Give a new free name to each suffixed row of a duplicate name.

    The row that got a suffix has a ``Non-unique Gene name`` (the base name) that is
    different from its ``Gene`` name. The other row owns the name and keeps it.
    The new name is ``<base>_N`` with the smallest free ``N >= 2``, as in the
    fixed genomes-catalogue-pipeline.

    :param rows: The CSV rows.
    :param duplicates: Output of :func:`find_duplicates` for the CSV names.
    :return: {row index: new name}.
    :raises PangenomeFixError: If a duplicate name has no owner row or more than one.
    """
    taken = {row[0] for row in rows}
    new_names = {}
    for name in sorted(duplicates):
        suffixed = [i for i in duplicates[name] if rows[i][1] not in ("", name)]
        owners = [i for i in duplicates[name] if i not in suffixed]
        if len(owners) != 1:
            raise PangenomeFixError(f"Name {name}: expected 1 row that owns the name, found {len(owners)}")
        for index in suffixed:
            base = rows[index][1]
            suffix = 2
            while f"{base}_{suffix}" in taken:
                suffix += 1
            new_names[index] = f"{base}_{suffix}"
            taken.add(new_names[index])
    return new_names


def extract_cds(gff: Path, wanted: dict[str, set[int]]) -> dict[str, str]:
    """Extract the CDS sequences of some locus tags from a Prokka GFF with a ``##FASTA`` section.

    Only the CDS with a length in the allowed set are extracted, because only those can
    be equal to a pan-genome record.

    :param gff: The GFF file. It can be compressed (``.gz``).
    :param wanted: {locus tag: allowed CDS lengths}.
    :return: {locus tag: CDS sequence (upper case, coding strand)}.
    :raises PangenomeFixError: If the GFF has candidate CDS but no ``##FASTA`` section.
    """
    features = {}  # locus tag -> (contig, start, end, strand)
    contigs = {}
    contig, chunks, in_fasta = None, [], False
    with open_text(gff) as handle:
        for line in handle:
            if not in_fasta:
                if line.startswith("##FASTA"):
                    in_fasta = True
                    needed_contigs = {feature[0] for feature in features.values()}
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) != 9 or parts[2] != "CDS":
                    continue
                attributes = dict(item.split("=", 1) for item in parts[8].split(";") if "=" in item)
                locus_tag = attributes.get("ID", "")
                start, end = int(parts[3]), int(parts[4])
                if locus_tag in wanted and end - start + 1 in wanted[locus_tag]:
                    features[locus_tag] = (parts[0], start, end, parts[6])
            elif line.startswith(">"):
                if contig in needed_contigs:
                    contigs[contig] = "".join(chunks)
                contig, chunks = line[1:].split()[0], []
            elif contig in needed_contigs:
                chunks.append(line.strip())
    if contig is not None and contig in needed_contigs:
        contigs[contig] = "".join(chunks)
    if features and not in_fasta:
        raise PangenomeFixError(f"{gff} has no ##FASTA section")

    sequences = {}
    for locus_tag, (contig, start, end, strand) in features.items():
        sequence = contigs[contig][start - 1 : end].upper()
        if strand == "-":
            sequence = sequence.translate(COMPLEMENT)[::-1]
        sequences[locus_tag] = sequence
    return sequences


def pair_duplicate_records(
    rows: list[list[str]],
    records: list[tuple[str, str]],
    duplicates: dict[str, list[int]],
    gff_dir: Path,
) -> dict[int, int]:
    """Pair each duplicate FASTA record with its CSV row.

    Each record is the sequence of the mmseqs representative of its cluster, and the
    representative is a member of the cluster. So a record belongs to the row that has a
    member CDS equal to the record sequence.

    :param rows: The CSV rows.
    :param records: The FASTA records, in file order.
    :param duplicates: Output of :func:`find_duplicates` for the CSV names.
    :param gff_dir: The folder with the member GFF files (``<genome>.gff.gz``).
    :return: {FASTA record index: CSV row index} for each duplicate record.
    :raises PangenomeFixError: If a record does not match exactly one row, or a GFF is missing.
    """
    record_positions = find_duplicates([name for name, _ in records])
    wanted = defaultdict(set)  # locus tag -> allowed lengths
    for name, row_indexes in duplicates.items():
        lengths = {len(records[i][1]) for i in record_positions.get(name, [])}
        for row_index in row_indexes:
            for cell in rows[row_index][CSV_META_COLUMNS:]:
                for locus_tag in split_cell(cell):
                    wanted[locus_tag] |= lengths

    by_genome = defaultdict(dict)
    for locus_tag, lengths in wanted.items():
        by_genome[genome_of(locus_tag)][locus_tag] = lengths

    member_sequences = {}
    for genome in sorted(by_genome):
        gff = gff_dir / f"{genome}.gff.gz"
        if not gff.exists():
            gff = gff_dir / f"{genome}.gff"
        if not gff.exists():
            raise PangenomeFixError(f"Member GFF not found for {genome} in {gff_dir}")
        member_sequences.update(extract_cds(gff, by_genome[genome]))

    pairs = {}
    for name, row_indexes in duplicates.items():
        record_indexes = record_positions.get(name, [])
        if len(record_indexes) != len(row_indexes):
            raise PangenomeFixError(f"Name {name}: {len(row_indexes)} CSV rows but {len(record_indexes)} FASTA records")
        row_sequences = {
            row_index: {
                member_sequences[tag]
                for cell in rows[row_index][CSV_META_COLUMNS:]
                for tag in split_cell(cell)
                if tag in member_sequences
            }
            for row_index in row_indexes
        }
        for record_index in record_indexes:
            sequence = records[record_index][1].upper()
            matches = [row_index for row_index in row_indexes if sequence in row_sequences[row_index]]
            if len(matches) != 1:
                raise PangenomeFixError(
                    f"Name {name}: FASTA record {record_index} matches {len(matches)} CSV rows, expected 1"
                )
            pairs[record_index] = matches[0]
        if sorted(pairs[i] for i in record_indexes) != sorted(row_indexes):
            raise PangenomeFixError(f"Name {name}: two FASTA records match the same CSV row")
    return pairs


def write_csv(path: Path, header: list[str], rows: list[list[str]]) -> None:
    """Write the CSV with paralogs joined by ``;``.

    :param path: The output file (compressed).
    :param header: The CSV header.
    :param rows: The CSV rows.
    """
    with open_gz_writer(path) as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(header)
        for row in rows:
            cells = [";".join(split_cell(cell)) for cell in row[CSV_META_COLUMNS:]]
            writer.writerow(row[:CSV_META_COLUMNS] + cells)


def write_rtab(path: Path, header: list[str], rows: list[list[str]]) -> None:
    """Write an Rtab built from the CSV, sorted by gene name as in the catalogue pipeline.

    :param path: The output file (compressed).
    :param header: The CSV header.
    :param rows: The CSV rows.
    """
    with open_gz_writer(path) as handle:
        handle.write("\t".join(["Gene"] + header[CSV_META_COLUMNS:]) + "\n")
        for row in sorted(rows, key=lambda row: row[0]):
            values = ["1" if split_cell(cell) else "0" for cell in row[CSV_META_COLUMNS:]]
            handle.write("\t".join([row[0]] + values) + "\n")


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    """Write the FASTA records, one sequence line per record.

    :param path: The output file (compressed).
    :param records: The (name, sequence) records.
    """
    with open_gz_writer(path) as handle:
        for name, sequence in records:
            handle.write(f">{name}\n{sequence}\n")


def fix_pangenome(pangenome_dir: Path, gff_dir: Path, species: str, outdir: Path) -> dict[str, str]:
    """Check the pangenome files of one species and write corrected files if needed.

    :param pangenome_dir: The ``pan-genome`` folder of the species.
    :param gff_dir: The folder with the member GFF files.
    :param species: The species representative accession (used in the report).
    :param outdir: The output folder.
    :return: The report line as {column: value}.
    :raises PangenomeFixError: If the files are not consistent and cannot be corrected.
    """
    header, rows = read_csv(pangenome_dir / CSV_NAME)
    rtab_rows = count_rtab_rows(pangenome_dir / RTAB_NAME)
    records = read_fasta(pangenome_dir / FASTA_NAME)
    names = [row[0] for row in rows]
    duplicates = find_duplicates(names)
    fasta_duplicates = find_duplicates([name for name, _ in records])

    report = {
        "species": species,
        "status": "ok",
        "csv_rows": str(len(rows)),
        "rtab_rows": str(rtab_rows),
        "fasta_records": str(len(records)),
        "duplicate_names": str(len(duplicates)),
        "renamed_rows": "0",
    }
    if len(rows) == rtab_rows == len(records) and not duplicates and not fasta_duplicates:
        return report
    if not duplicates or Counter(names) != Counter(name for name, _ in records):
        raise PangenomeFixError(
            f"{species}: files are not consistent, but not because of duplicate names "
            f"(CSV {len(rows)}, Rtab {rtab_rows}, FASTA {len(records)})"
        )

    new_names = new_names_for_duplicates(rows, duplicates)
    pairs = pair_duplicate_records(rows, records, duplicates, gff_dir)

    renames = []
    for index, new_name in sorted(new_names.items()):
        renames.append((rows[index][0], new_name, rows[index][1]))
        rows[index][0] = new_name
    fixed_records = [
        (rows[pairs[i]][0], sequence) if i in pairs else (name, sequence) for i, (name, sequence) in enumerate(records)
    ]
    if len({name for name, _ in fixed_records}) != len(rows):
        raise PangenomeFixError(f"{species}: corrected names are not unique")

    outdir.mkdir(parents=True, exist_ok=True)
    write_csv(outdir / CSV_NAME, header, rows)
    write_rtab(outdir / RTAB_NAME, header, rows)
    write_fasta(outdir / FASTA_NAME, fixed_records)
    with open(outdir / "renamed_genes.tsv", "w") as handle:
        handle.write("old_name\tnew_name\tnon_unique_gene_name\n")
        for old_name, new_name, base in renames:
            handle.write(f"{old_name}\t{new_name}\t{base}\n")

    report["status"] = "fixed"
    report["renamed_rows"] = str(len(renames))
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--pangenome_dir", type=Path, required=True, help="The pan-genome folder of the species")
    parser.add_argument("--gff_dir", type=Path, required=True, help="The folder with the member GFF files")
    parser.add_argument("--species", required=True, help="The species representative accession")
    parser.add_argument("--outdir", type=Path, default=Path("fixed"), help="The folder for the corrected files")
    parser.add_argument("--report", type=Path, required=True, help="The report TSV (one line per species)")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    try:
        report = fix_pangenome(args.pangenome_dir, args.gff_dir, args.species, args.outdir)
    except PangenomeFixError as error:
        logging.error(error)
        sys.exit(1)

    with open(args.report, "w") as handle:
        handle.write("\t".join(report) + "\n")
        handle.write("\t".join(report.values()) + "\n")
    logging.info("%s: %s (%s renamed rows)", args.species, report["status"], report["renamed_rows"])


if __name__ == "__main__":
    main()
