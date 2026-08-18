#!/usr/bin/env python3
"""Build a PolyMarker marker list from KASPberry SNPs."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a PolyMarker marker list."
    )
    parser.add_argument("--snps", required=True, type=Path)
    parser.add_argument("--snp-status", required=True, type=Path)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--aliases", required=True, type=Path)
    parser.add_argument("--flank", required=True, type=int)
    parser.add_argument("--out", required=True, type=Path)

    return parser.parse_args()


def read_reference(path: Path) -> str:
    """Read the single-sequence regional FASTA."""
    with path.open(encoding="utf-8") as handle:
        return "".join(
            line.strip()
            for line in handle
            if not line.startswith(">")
        ).upper()


def read_target(aliases_path: Path) -> tuple[str, str]:
    """Return the genotype and PolyMarker target chromosome."""
    with aliases_path.open(newline="", encoding="utf-8") as handle:
        source = next(
            row
            for row in csv.DictReader(handle, delimiter="\t")
            if row["is_source_seq"] == "true"
        )

    return source["genotype"], source["alias"][:2]


def read_passed_snp_ids(path: Path) -> set[str]:
    """Return SNP IDs that passed diagnostic filtering."""
    with path.open(newline="", encoding="utf-8") as handle:
        return {
            row["snp_id"]
            for row in csv.DictReader(handle, delimiter="\t")
            if row["status"] == "PASS"
        }


def read_snps(
    snps_path: Path,
    genotype: str,
    passed_snp_ids: set[str],
) -> list[tuple[str, int, str, str]]:
    """Return SNPs and alleles for one genotype."""
    grouped = {}

    with snps_path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["snp_id"] not in passed_snp_ids:
                continue

            grouped.setdefault(row["snp_id"], []).append(row)

    snps = []

    for snp_id, rows in grouped.items():
        genotype_row = next(
            row
            for row in rows
            if row["sample"] == genotype
        )

        allele = genotype_row["nt"].upper()

        other_allele = next(
            row["nt"].upper()
            for row in rows
            if row["nt"].upper() != allele
        )

        snps.append(
            (
                snp_id,
                int(genotype_row["pos_in_zone"]),
                allele,
                other_allele,
            )
        )

    return snps


def main() -> None:
    args = parse_args()

    sequence = read_reference(args.reference)
    genotype, target_chr = read_target(args.aliases)
    passed_snp_ids = read_passed_snp_ids(args.snp_status)
    snps = read_snps(args.snps, genotype, passed_snp_ids)

    args.out.parent.mkdir(parents=True, exist_ok=True)

    with args.out.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.writer(
            handle,
            lineterminator="\n",
        )

        for marker_id, pos1, allele, other_allele in snps:
            pos0 = pos1 - 1

            start = max(0, pos0 - args.flank)
            end = min(len(sequence), pos0 + args.flank + 1)

            left = sequence[start:pos0]
            right = sequence[pos0 + 1:end]

            marker_sequence = (
                f"{left}[{allele}/{other_allele}]{right}"
            )

            writer.writerow([
                marker_id,
                target_chr,
                marker_sequence,
            ])

    print(f"{genotype}: wrote {len(snps)} PolyMarker markers")


if __name__ == "__main__":
    main()
