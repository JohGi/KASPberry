#!/usr/bin/env python3
# Author: Johanna Girodolle

"""Write block coordinates in region and source-sequence coordinate systems."""

from __future__ import annotations

import argparse
import csv
import logging
from pathlib import Path

import polars as pl
from attrs import define, field

LOGGER = logging.getLogger(__name__)


@define(frozen=True)
class GenotypeOffset:
    """Store the source-sequence start offset for one genotype."""

    genotype: str
    region_start_in_source_seq: int = 1


@define(frozen=True)
class BlockRecord:
    """Store one block coordinate record for one genotype."""

    block_id: str
    genotype: str
    block_start_in_region: int
    block_end_in_region: int
    block_start_in_source_seq: int
    block_end_in_source_seq: int


@define
class BlockCoordinateWriter:
    """Build block coordinate records from a GFF and a genotypes TSV."""

    gff_path: Path
    genotypes_tsv_path: Path
    output_path: Path
    genotype_offsets: dict[str, GenotypeOffset] = field(factory=dict)

    def run(self) -> None:
        """Run the full block coordinate export workflow."""
        self.genotype_offsets = read_genotype_offsets(self.genotypes_tsv_path)
        block_records = read_block_records(self.gff_path, self.genotype_offsets)
        dataframe = build_block_dataframe(block_records)
        write_dataframe(dataframe, self.output_path)


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Write block coordinates from filtered block GFF coordinates in the region "
            "to both region and source-sequence coordinate systems."
        )
    )
    parser.add_argument(
        "--gff",
        required=True,
        type=Path,
        help="Filtered block GFF file.",
    )
    parser.add_argument(
        "--genotypes-tsv",
        required=True,
        type=Path,
        help=(
            "Genotypes TSV. Column 2 must contain genotype names. "
            "Column 3 is optional and, if present, is interpreted as "
            "region_start_in_source_seq. Missing or empty offsets default to 1."
        ),
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Output TSV path.",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging level.",
    )
    return parser.parse_args()


def setup_logging(log_level: str) -> None:
    """Configure logging."""
    logging.basicConfig(
        level=getattr(logging, log_level),
        format="%(levelname)s | %(message)s",
    )


def compute_projected_position(start_position: int, relative_position: int) -> int:
    """Project a 1-based relative position onto a 1-based coordinate system."""
    return start_position + relative_position - 1


def natural_sort_key(value: str) -> tuple[int, str]:
    """Sort numerically when possible, otherwise lexicographically."""
    if value.isdigit():
        return int(value), value
    return 10**18, value


def extract_block_id(attributes: str, gff_path: Path, line_number: int) -> str:
    """Extract the block ID from the GFF attributes column."""
    for attribute in attributes.split(";"):
        if attribute.startswith("ID="):
            return attribute.removeprefix("ID=")

    raise ValueError(
        f"Could not find ID attribute in column 9 of {gff_path} at line {line_number}: "
        f"{attributes}"
    )


def read_genotype_offsets(genotypes_tsv_path: Path) -> dict[str, GenotypeOffset]:
    """Read per-genotype source-sequence offsets from the genotypes TSV."""
    genotype_offsets: dict[str, GenotypeOffset] = {}

    with genotypes_tsv_path.open("r", encoding="utf-8", newline="") as handle:
        for line_number, row in enumerate(
            csv.DictReader(handle, delimiter="\t"),
            start=2,
        ):
            genotype = (row.get("genotype") or "").strip()

            if not genotype:
                raise ValueError(
                    f"Empty genotype name in {genotypes_tsv_path} "
                    f"at line {line_number}"
                )

            value = (row.get("region_start") or "").strip()
            region_start_in_source_seq = int(value) if value else 1

            genotype_offsets[genotype] = GenotypeOffset(
                genotype=genotype,
                region_start_in_source_seq=region_start_in_source_seq,
            )

    LOGGER.info(
        "Read %d genotype offsets from %s",
        len(genotype_offsets),
        genotypes_tsv_path,
    )

    return genotype_offsets


def read_block_records(
    gff_path: Path,
    genotype_offsets: dict[str, GenotypeOffset],
) -> list[BlockRecord]:
    """Read block records from a GFF and project them to source coordinates."""
    block_records: list[BlockRecord] = []

    with gff_path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(
                    f"Expected 9 tab-separated columns in GFF {gff_path} at line {line_number}, "
                    f"got {len(fields)}: {line.rstrip()}"
                )

            genotype = fields[0]
            block_start_in_region = int(fields[3])
            block_end_in_region = int(fields[4])
            block_id = extract_block_id(fields[8], gff_path, line_number)

            region_start_in_source_seq = genotype_offsets.get(genotype, GenotypeOffset(genotype)).region_start_in_source_seq

            block_start_in_source_seq = compute_projected_position(
                region_start_in_source_seq,
                block_start_in_region,
            )
            block_end_in_source_seq = compute_projected_position(
                region_start_in_source_seq,
                block_end_in_region,
            )

            block_records.append(
                BlockRecord(
                    block_id=block_id,
                    genotype=genotype,
                    block_start_in_region=block_start_in_region,
                    block_end_in_region=block_end_in_region,
                    block_start_in_source_seq=block_start_in_source_seq,
                    block_end_in_source_seq=block_end_in_source_seq,
                )
            )

    LOGGER.info("Read %d block records from %s", len(block_records), gff_path)
    return block_records


def build_block_dataframe(block_records: list[BlockRecord]) -> pl.DataFrame:
    """Build the output dataframe."""
    rows = [
        {
            "block_id": record.block_id,
            "genotype": record.genotype,
            "block_start_in_region": record.block_start_in_region,
            "block_end_in_region": record.block_end_in_region,
            "block_start_in_source_seq": record.block_start_in_source_seq,
            "block_end_in_source_seq": record.block_end_in_source_seq,
        }
        for record in block_records
    ]

    dataframe = pl.DataFrame(rows)

    sorted_rows = sorted(
        dataframe.iter_rows(named=True),
        key=lambda row: (natural_sort_key(str(row["block_id"])), str(row["genotype"])),
    )
    return pl.DataFrame(sorted_rows)


def write_dataframe(dataframe: pl.DataFrame, output_path: Path) -> None:
    """Write a dataframe as a TSV file."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    dataframe.write_csv(output_path, separator="\t")
    LOGGER.info("Wrote %s", output_path)


def main() -> None:
    """Run the block coordinate export script."""
    args = parse_args()
    setup_logging(args.log_level)

    writer = BlockCoordinateWriter(
        gff_path=args.gff,
        genotypes_tsv_path=args.genotypes_tsv,
        output_path=args.output,
    )
    writer.run()


if __name__ == "__main__":
    main()
