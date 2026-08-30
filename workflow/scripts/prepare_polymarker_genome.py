#!/usr/bin/env python3
"""Prepare a PolyMarker-compatible FASTA for one genotype."""

from __future__ import annotations

import argparse
import csv
import string
from pathlib import Path


GROUP_CODE_CHARS = (
    string.digits
    + string.ascii_uppercase
    + string.ascii_lowercase
    + "-_.:*#"
)

SUBGENOME_CODE_CHARS = (
    string.digits
    + string.ascii_uppercase
    + string.ascii_lowercase
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prepare a PolyMarker-compatible genome FASTA."
    )
    parser.add_argument("--genotype", required=True)
    parser.add_argument("--source-seq", required=True)
    parser.add_argument("--genome", required=True, type=Path)
    parser.add_argument("--chromosomes", default="")
    parser.add_argument("--out-fasta", required=True, type=Path)
    parser.add_argument("--out-aliases", required=True, type=Path)
    return parser.parse_args()


def read_fasta_ids(path: Path) -> list[str]:
    """Return FASTA record IDs in input order."""
    ids = []

    with path.open("rb") as handle:
        for line in handle:
            if line.startswith(b">"):
                ids.append(
                    line[1:]
                    .split(None, 1)[0]
                    .decode("utf-8")
                )

    return ids


def assign_codes(
    labels: list[str],
    code_chars: str,
) -> dict[str, str]:
    """Assign one-character codes, preserving simple labels when possible."""
    labels = sorted(set(labels))
    codes = {}
    used = set()

    for label in labels:
        if len(label) == 1 and label in code_chars:
            codes[label] = label
            used.add(label)

    available = iter(
        char
        for char in code_chars
        if char not in used
    )

    for label in labels:
        if label not in codes:
            codes[label] = next(available)

    return codes


def build_aliases(
    genome_ids: list[str],
    chromosomes_path: Path | None,
    genotype: str,
) -> tuple[dict[str, str], dict[str, tuple[str, str]]]:
    """Build PolyMarker aliases for every genome FASTA record."""

    # No chromosome information: every FASTA record is treated
    # as an independent chromosome group.
    if chromosomes_path is None:
        aliases = {
            seq_id: f"{group}A__001"
            for seq_id, group in zip(genome_ids, GROUP_CODE_CHARS)
        }
        metadata = {
            seq_id: ("", "")
            for seq_id in genome_ids
        }
        return aliases, metadata

    with chromosomes_path.open(
        newline="",
        encoding="utf-8",
    ) as handle:
        rows = [
            row
            for row in csv.DictReader(handle, delimiter="\t")
            if row["genotype"] == genotype
        ]

    group_codes = assign_codes(
        [
            row["homoeologous_group"]
            for row in rows
        ],
        GROUP_CODE_CHARS,
    )

    subgenome_codes = assign_codes(
        [
            row["subgenome"]
            for row in rows
        ],
        SUBGENOME_CODE_CHARS,
    )

    aliases = {}
    metadata = {}
    counters = {}

    for row in sorted(
        rows,
        key=lambda row: (
            row["homoeologous_group"],
            row["subgenome"],
            row["seq_id"],
        ),
    ):
        seq_id = row["seq_id"]
        group = row["homoeologous_group"]
        subgenome = row["subgenome"]

        pair = (group, subgenome)
        counters[pair] = counters.get(pair, 0) + 1

        aliases[seq_id] = (
            f"{group_codes[group]}"
            f"{subgenome_codes[subgenome]}"
            f"__{counters[pair]:03d}"
        )

        metadata[seq_id] = (group, subgenome)

    # All records absent from chromosomes.tsv are retained
    # in one artificial, unassigned PolyMarker chromosome.
    unassigned = [
        seq_id
        for seq_id in genome_ids
        if seq_id not in aliases
    ]

    if unassigned:
        used_groups = set(group_codes.values())
        unassigned_group = next(
            char
            for char in GROUP_CODE_CHARS
            if char not in used_groups
        )

        for i, seq_id in enumerate(unassigned, start=1):
            aliases[seq_id] = f"{unassigned_group}0__{i:06d}"
            metadata[seq_id] = ("", "")

    return aliases, metadata


def write_polymarker_fasta(
    genome_path: Path,
    output_path: Path,
    aliases: dict[str, str],
) -> None:
    """Write the complete genome with PolyMarker-compatible record IDs."""
    with (
        genome_path.open("rb", buffering=16 * 1024 * 1024) as source,
        output_path.open("wb", buffering=16 * 1024 * 1024) as output,
    ):
        for line in source:
            if line.startswith(b">"):
                seq_id = (
                    line[1:]
                    .split(None, 1)[0]
                    .decode("utf-8")
                )
                output.write(
                    f">{aliases[seq_id]}\n".encode("ascii")
                )
            else:
                output.write(line)


def write_alias_table(
    output_path: Path,
    genotype: str,
    source_seq: str,
    aliases: dict[str, str],
    metadata: dict[str, tuple[str, str]],
) -> None:
    """Write original-to-PolyMarker sequence ID mapping."""
    with output_path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.writer(
            handle,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writerow([
            "genotype",
            "original_seq_id",
            "alias",
            "homoeologous_group",
            "subgenome",
            "is_source_seq",
        ])

        for seq_id in sorted(aliases):
            group, subgenome = metadata[seq_id]

            writer.writerow([
                genotype,
                seq_id,
                aliases[seq_id],
                group,
                subgenome,
                str(seq_id == source_seq).lower(),
            ])


def main() -> None:
    args = parse_args()

    chromosomes_path = (
        Path(args.chromosomes)
        if args.chromosomes
        else None
    )

    genome_ids = read_fasta_ids(args.genome)

    aliases, metadata = build_aliases(
        genome_ids=genome_ids,
        chromosomes_path=chromosomes_path,
        genotype=args.genotype,
    )

    args.out_fasta.parent.mkdir(parents=True, exist_ok=True)
    args.out_aliases.parent.mkdir(parents=True, exist_ok=True)

    write_polymarker_fasta(
        args.genome,
        args.out_fasta,
        aliases,
    )

    write_alias_table(
        args.out_aliases,
        args.genotype,
        args.source_seq,
        aliases,
        metadata,
    )

    print(
        f"{args.genotype}: prepared {len(aliases)} genome sequences "
        "for PolyMarker and MFEprimer"
    )


if __name__ == "__main__":
    main()
