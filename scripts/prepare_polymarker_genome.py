#!/usr/bin/env python3
"""Prepare a PolyMarker-compatible FASTA for one genotype."""

from __future__ import annotations

import argparse
import csv
import string
from pathlib import Path


CODE_CHARS = string.digits + string.ascii_uppercase + string.ascii_lowercase


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


def assign_codes(labels: list[str]) -> dict[str, str]:
    """Assign one-character codes, preserving simple labels when possible."""
    labels = sorted(set(labels))
    codes = {}
    used = set()

    for label in labels:
        if len(label) == 1 and label in CODE_CHARS:
            codes[label] = label
            used.add(label)

    available = iter(
        char for char in CODE_CHARS
        if char not in used
    )

    for label in labels:
        if label not in codes:
            try:
                codes[label] = next(available)
            except StopIteration as error:
                raise ValueError(
                    "Too many labels to encode in PolyMarker sequence IDs."
                ) from error

    return codes


def build_aliases(
    chromosomes_path: Path | None,
    genotype: str,
    source_seq: str,
) -> tuple[dict[str, str], dict[str, tuple[str, str]]]:
    """Build PolyMarker aliases for the sequences to retain."""

    if chromosomes_path is None:
        return (
            {source_seq: "0A__001"},
            {source_seq: ("", "A")},
        )

    with chromosomes_path.open(
        newline="",
        encoding="utf-8",
    ) as handle:
        rows = [
            row
            for row in csv.DictReader(handle, delimiter="\t")
            if row["genotype"] == genotype
        ]

    group_codes = assign_codes([
        row["homoeologous_group"]
        for row in rows
    ])

    subgenome_codes = assign_codes([
        row["subgenome"]
        for row in rows
    ])

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

    return aliases, metadata


def write_polymarker_fasta(
    genome_path: Path,
    output_path: Path,
    aliases: dict[str, str],
) -> None:
    """Stream the genome and write only sequences used by PolyMarker."""

    seen = set()
    keep = False

    with (
        genome_path.open(
            "rb",
            buffering=16 * 1024 * 1024,
        ) as source,
        output_path.open(
            "wb",
            buffering=16 * 1024 * 1024,
        ) as output,
    ):
        for line in source:
            if line.startswith(b">"):
                seq_id = (
                    line[1:]
                    .split(None, 1)[0]
                    .decode("utf-8")
                )

                keep = seq_id in aliases

                if keep:
                    seen.add(seq_id)
                    output.write(
                        f">{aliases[seq_id]}\n".encode("ascii")
                    )

            elif keep:
                output.write(line)

    missing = set(aliases) - seen

    if missing:
        raise ValueError(
            "Sequences expected by PolyMarker are missing from "
            f"{genome_path}: "
            + ", ".join(sorted(missing))
        )


def write_alias_table(
    output_path: Path,
    genotype: str,
    source_seq: str,
    aliases: dict[str, str],
    metadata: dict[str, tuple[str, str]],
) -> None:
    """Write the mapping between original and PolyMarker sequence IDs."""

    with output_path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")

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

    aliases, metadata = build_aliases(
        chromosomes_path=chromosomes_path,
        genotype=args.genotype,
        source_seq=args.source_seq,
    )

    args.out_fasta.parent.mkdir(
        parents=True,
        exist_ok=True,
    )
    args.out_aliases.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    write_polymarker_fasta(
        genome_path=args.genome,
        output_path=args.out_fasta,
        aliases=aliases,
    )

    write_alias_table(
        output_path=args.out_aliases,
        genotype=args.genotype,
        source_seq=args.source_seq,
        aliases=aliases,
        metadata=metadata,
    )

    print(
        f"{args.genotype}: retained {len(aliases)} sequences "
        "for PolyMarker"
    )


if __name__ == "__main__":
    main()
