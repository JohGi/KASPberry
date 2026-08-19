#!/usr/bin/env python3
# Author: Johanna Girodolle

"""Write summary statistics for filtered blocks and SNPs."""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

import polars as pl
from attrs import define

LOGGER = logging.getLogger(__name__)


@define(frozen=True)
class GenotypeStats:
    """Store summary statistics for one genotype."""

    genotype: str
    region_length_bp: int
    n_blocks_present: int
    cumulated_block_bp: int
    covered_pct_of_region: float

    def to_dict(self) -> dict[str, int | float]:
        """Convert the genotype statistics to a serializable dictionary."""
        return {
            "region_length_bp": self.region_length_bp,
            "n_blocks_present": self.n_blocks_present,
            "cumulated_block_bp": self.cumulated_block_bp,
            "covered_pct_of_region": round(self.covered_pct_of_region, 2),
        }


@define(frozen=True)
class GlobalStats:
    """Store global summary statistics."""

    n_blocks_kept: int
    min_block_len_bp: int
    max_block_len_bp: int
    mean_block_len_bp: float
    n_snps_kept: int

    def to_dict(self) -> dict[str, int | float]:
        """Convert the global statistics to a serializable dictionary."""
        return {
            "n_blocks_kept": self.n_blocks_kept,
            "min_block_len_bp": self.min_block_len_bp,
            "max_block_len_bp": self.max_block_len_bp,
            "mean_block_len_bp": round(self.mean_block_len_bp, 2),
            "n_snps_kept": self.n_snps_kept,
        }


@define(frozen=True)
class SummaryStats:
    """Store all summary statistics."""

    global_stats: GlobalStats
    genotype_stats: dict[str, GenotypeStats]

    def to_dict(self) -> dict[str, dict]:
        """Convert summary statistics to a serializable dictionary."""
        return {
            "global": self.global_stats.to_dict(),
            "samples": {
                genotype: stats.to_dict()
                for genotype, stats in sorted(self.genotype_stats.items())
            },
        }


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Write summary statistics for filtered blocks and SNPs."
    )
    parser.add_argument(
        "--block-coords",
        required=True,
        help="TSV file containing filtered block coordinates.",
    )
    parser.add_argument(
        "--snp-positions",
        required=True,
        help="TSV file containing SNP positions in wide format.",
    )
    parser.add_argument(
        "--clean-fastas",
        required=True,
        nargs="+",
        help="Clean FASTA files, one per genotype.",
    )
    parser.add_argument(
        "--json-output",
        required=True,
        help="Output JSON file with machine-readable summary statistics.",
    )
    parser.add_argument(
        "--txt-output",
        required=True,
        help="Output TXT file with human-readable summary statistics.",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Logging level.",
    )
    return parser.parse_args()


def configure_logging(log_level: str) -> None:
    """Configure application logging."""
    logging.basicConfig(
        level=getattr(logging, log_level),
        format="[%(levelname)s] %(message)s",
    )


def genotype_name_from_fasta_path(fasta_path: Path) -> str:
    """Infer the genotype name from a FASTA filename stem."""
    return fasta_path.stem


def read_single_fasta_length(fasta_path: Path) -> int:
    """Return the sequence length of a single-record FASTA file."""
    sequence_parts: list[str] = []

    with fasta_path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith(">"):
                continue
            sequence_parts.append(line)

    return len("".join(sequence_parts))


def read_fasta_lengths(fasta_paths: list[Path]) -> dict[str, int]:
    """Read clean FASTA lengths for all genotypes."""
    return {
        genotype_name_from_fasta_path(fasta_path): read_single_fasta_length(fasta_path)
        for fasta_path in fasta_paths
    }


def read_block_coordinates(block_coords_path: Path) -> pl.DataFrame:
    """Read block coordinates and compute block lengths."""
    return pl.read_csv(block_coords_path, separator="\t").with_columns(
        (
            pl.col("block_end_in_region") - pl.col("block_start_in_region") + 1
        ).alias("block_len_bp")
    )


def read_snp_positions(snp_positions_path: Path) -> pl.DataFrame:
    """Read the SNP positions table."""
    return pl.read_csv(snp_positions_path, separator="\t")


def compute_global_stats(
    block_coords_df: pl.DataFrame,
    snp_positions_df: pl.DataFrame,
) -> GlobalStats:
    """Compute global summary statistics."""
    stats_row = block_coords_df.select(
        [
            pl.col("block_id").n_unique().alias("n_blocks_kept"),
            pl.col("block_len_bp").min().alias("min_block_len_bp"),
            pl.col("block_len_bp").max().alias("max_block_len_bp"),
            pl.col("block_len_bp").mean().alias("mean_block_len_bp"),
        ]
    ).to_dicts()[0]

    return GlobalStats(
        n_blocks_kept=int(stats_row["n_blocks_kept"]),
        min_block_len_bp=int(stats_row["min_block_len_bp"]),
        max_block_len_bp=int(stats_row["max_block_len_bp"]),
        mean_block_len_bp=float(stats_row["mean_block_len_bp"]),
        n_snps_kept=int(snp_positions_df.height),
    )


def compute_genotype_stats(
    block_coords_df: pl.DataFrame,
    fasta_lengths: dict[str, int],
) -> dict[str, GenotypeStats]:
    """Compute per-genotype summary statistics."""
    genotype_rows = block_coords_df.group_by("genotype").agg(
        [
            pl.col("block_id").n_unique().alias("n_blocks_present"),
            pl.col("block_len_bp").sum().alias("cumulated_block_bp"),
        ]
    )

    stats_by_genotype: dict[str, GenotypeStats] = {}

    for row in genotype_rows.iter_rows(named=True):
        genotype = str(row["genotype"])
        region_length_bp = fasta_lengths[genotype]
        cumulated_block_bp = int(row["cumulated_block_bp"])
        covered_pct_of_region = (cumulated_block_bp / region_length_bp) * 100

        stats_by_genotype[genotype] = GenotypeStats(
            genotype=genotype,
            region_length_bp=region_length_bp,
            n_blocks_present=int(row["n_blocks_present"]),
            cumulated_block_bp=cumulated_block_bp,
            covered_pct_of_region=covered_pct_of_region,
        )

    return stats_by_genotype


def build_summary_stats(
    block_coords_path: Path,
    snp_positions_path: Path,
    clean_fasta_paths: list[Path],
) -> SummaryStats:
    """Build summary statistics from workflow outputs."""
    block_coords_df = read_block_coordinates(block_coords_path)
    snp_positions_df = read_snp_positions(snp_positions_path)
    fasta_lengths = read_fasta_lengths(clean_fasta_paths)

    return SummaryStats(
        global_stats=compute_global_stats(block_coords_df, snp_positions_df),
        genotype_stats=compute_genotype_stats(block_coords_df, fasta_lengths),
    )


def write_summary_json(summary_stats: SummaryStats, output_path: Path) -> None:
    """Write summary statistics to JSON."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as handle:
        json.dump(summary_stats.to_dict(), handle, indent=2)
        handle.write("\n")


def build_summary_text(summary_stats: SummaryStats) -> str:
    """Build a human-readable summary text."""
    lines = [
        "KASPberry run summary",
        "=====================",
        "",
        "Global",
        "------",
        f"Kept blocks: {summary_stats.global_stats.n_blocks_kept}",
        f"Smallest block length (bp): {summary_stats.global_stats.min_block_len_bp}",
        f"Largest block length (bp): {summary_stats.global_stats.max_block_len_bp}",
        f"Mean block length (bp): {summary_stats.global_stats.mean_block_len_bp:.2f}",
        f"Kept SNPs: {summary_stats.global_stats.n_snps_kept}",
        "",
        "Per-genotype",
        "------------",
    ]

    for genotype, stats in sorted(summary_stats.genotype_stats.items()):
        lines.extend(
            [
                genotype,
                f"  Region length (bp): {stats.region_length_bp}",
                f"  Blocks present: {stats.n_blocks_present}",
                f"  Cumulated block length (bp): {stats.cumulated_block_bp}",
                f"  Covered region (%): {stats.covered_pct_of_region:.2f}",
                "",
            ]
        )

    return "\n".join(lines).rstrip() + "\n"


def write_summary_text(summary_stats: SummaryStats, output_path: Path) -> None:
    """Write summary statistics to a human-readable text file."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(build_summary_text(summary_stats), encoding="utf-8")


def main() -> None:
    """Run the summary statistics writer."""
    args = parse_args()
    configure_logging(args.log_level)

    summary_stats = build_summary_stats(
        block_coords_path=Path(args.block_coords),
        snp_positions_path=Path(args.snp_positions),
        clean_fasta_paths=[Path(path) for path in args.clean_fastas],
    )
    write_summary_json(summary_stats, Path(args.json_output))
    write_summary_text(summary_stats, Path(args.txt_output))

    LOGGER.info("Summary JSON written to %s", args.json_output)
    LOGGER.info("Summary TXT written to %s", args.txt_output)


if __name__ == "__main__":
    main()
