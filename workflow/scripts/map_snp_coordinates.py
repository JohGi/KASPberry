#!/usr/bin/env python3
# Author: Johanna Girodolle

"""Map SNP coordinates from alignment coordinates to block, region, and source sequence coordinates."""

from __future__ import annotations

import argparse
import csv
import logging
from pathlib import Path

import polars as pl
from attrs import define, field
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

from ids import make_snp_id

LOGGER = logging.getLogger(__name__)


@define(frozen=True)
class GenotypeOffset:
    """Store the source-sequence start offset for one genotype."""

    genotype: str
    region_start_in_source_seq: int = 1


@define(frozen=True)
class VariantRecord:
    """Store one SNP record from the filtered VCF."""

    block_id: str
    aln_pos: int


@define(frozen=True)
class LongRow:
    """Store one long-format output row."""

    block_id: str
    aln_pos: int
    genotype: str
    nt: str
    pos_in_block: int | None
    block_start_in_region: int | None
    pos_in_region: int | None
    region_start_in_source_seq: int | None
    pos_in_source_seq: int | None

    @property
    def snp_id(self) -> str:
        """Return the canonical SNP identifier."""
        return make_snp_id(self.block_id, self.aln_pos)


@define(frozen=True)
class SnpProjection:
    """Store projected values for one genotype at one SNP position."""

    nt: str
    pos_in_block: int | None


@define
class AlignmentProjector:
    """Project selected alignment columns to ungapped positions for each genotype."""

    block_id: str
    alignment_path: Path
    genotype_order: list[str]
    projections_by_aln_pos: dict[int, dict[str, SnpProjection]] = field(factory=dict)

    def load(self, target_aln_positions: set[int]) -> None:
        """Load the alignment and precompute projections only for requested SNP columns."""
        if not target_aln_positions:
            self.projections_by_aln_pos = {}
            return

        alignment = AlignIO.read(str(self.alignment_path), "fasta")
        normalized_genotype_names = self.get_normalized_genotype_names(alignment)
        self.validate_alignment_genotypes(normalized_genotype_names)
        self.projections_by_aln_pos = self.build_projection_cache(
            alignment=alignment,
            normalized_genotype_names=normalized_genotype_names,
            target_aln_positions=target_aln_positions,
        )

    @staticmethod
    def normalize_alignment_genotype_name(sequence_id: str) -> str:
        """Extract the genotype name from an alignment record identifier."""
        return sequence_id.split(":", 1)[0]

    def get_normalized_genotype_names(self, alignment: MultipleSeqAlignment) -> list[str]:
        """Return normalized genotype names for all alignment records."""
        return [self.normalize_alignment_genotype_name(record.id) for record in alignment]

    def validate_alignment_genotypes(self, normalized_genotype_names: list[str]) -> None:
        """Validate alignment genotype names against the expected VCF genotype order."""
        observed_genotypes = set(normalized_genotype_names)
        expected_genotypes = set(self.genotype_order)

        if len(observed_genotypes) != len(normalized_genotype_names):
            raise ValueError(
                f"Duplicate normalized genotype names found in alignment for block {self.block_id}: "
                f"{normalized_genotype_names}"
            )

        if observed_genotypes != expected_genotypes:
            raise ValueError(
                f"Alignment genotypes do not match VCF genotypes for block {self.block_id}. "
                f"Expected={sorted(expected_genotypes)} Observed={sorted(observed_genotypes)}"
            )

    def build_projection_cache(
        self,
        alignment: MultipleSeqAlignment,
        normalized_genotype_names: list[str],
        target_aln_positions: set[int],
    ) -> dict[int, dict[str, SnpProjection]]:
        """Build per-genotype projections only for requested alignment positions."""
        alignment_length = alignment.get_alignment_length()
        max_target = max(target_aln_positions)
        ungapped_counters: dict[str, int] = {genotype: 0 for genotype in normalized_genotype_names}
        projections_by_aln_pos: dict[int, dict[str, SnpProjection]] = {}

        for aln_index in range(alignment_length):
            aln_pos = aln_index + 1

            for record, genotype in zip(alignment, normalized_genotype_names):
                nt = str(record.seq[aln_index]).upper()
                if nt != "-":
                    ungapped_counters[genotype] += 1

            if aln_pos not in target_aln_positions:
                if aln_pos >= max_target:
                    break
                continue

            projections_by_aln_pos[aln_pos] = {}

            for record, genotype in zip(alignment, normalized_genotype_names):
                nt = str(record.seq[aln_index]).upper()
                pos_in_block = ungapped_counters[genotype] if nt != "-" else None
                projections_by_aln_pos[aln_pos][genotype] = SnpProjection(
                    nt=nt,
                    pos_in_block=pos_in_block,
                )

            if aln_pos >= max_target and len(projections_by_aln_pos) == len(target_aln_positions):
                break

        missing_positions = sorted(target_aln_positions - set(projections_by_aln_pos))
        if missing_positions:
            raise ValueError(
                f"Alignment {self.alignment_path} does not contain requested SNP positions "
                f"for block {self.block_id}: {missing_positions}"
            )

        return projections_by_aln_pos

    def get_projection(self, aln_pos: int, genotype: str) -> SnpProjection:
        """Return the nucleotide and ungapped block position for one genotype at one SNP column."""
        return self.projections_by_aln_pos[aln_pos][genotype]


@define
class SnpPositionProjector:
    """Project SNP coordinates from alignments to region and source sequences."""

    vcf_path: Path
    block_coords_path: Path
    genotypes_tsv_path: Path
    align_dir: Path
    long_output_path: Path
    wide_output_path: Path
    genotype_order: list[str] = field(factory=list)
    variants_by_block: dict[str, list[VariantRecord]] = field(factory=dict)
    block_starts_in_region: dict[tuple[str, str], int] = field(factory=dict)
    genotype_offsets: dict[str, GenotypeOffset] = field(factory=dict)

    def run(self) -> None:
        """Run the full SNP projection workflow."""
        self.genotype_order, self.variants_by_block = read_vcf(self.vcf_path)
        self.block_starts_in_region = read_block_coords(self.block_coords_path)
        self.genotype_offsets = read_genotype_offsets(self.genotypes_tsv_path)
        long_rows = self.project_variants()
        long_df = build_long_dataframe(long_rows)
        wide_df = build_wide_dataframe(long_rows, self.genotype_order)
        write_dataframe(long_df, self.long_output_path)
        write_dataframe(wide_df, self.wide_output_path)

    def project_variants(self) -> list[LongRow]:
        """Project all variants to long-format rows."""
        long_rows: list[LongRow] = []

        for block_id in sorted(self.variants_by_block, key=natural_sort_key):
            alignment_path = self.align_dir / f"{block_id}.aln.fasta"
            block_variants = self.variants_by_block[block_id]
            target_aln_positions = {variant.aln_pos for variant in block_variants}

            LOGGER.info(
                "Projecting block %s from alignment %s using %d SNP columns",
                block_id,
                alignment_path,
                len(target_aln_positions),
            )

            projector = AlignmentProjector(
                block_id=block_id,
                alignment_path=alignment_path,
                genotype_order=self.genotype_order,
            )
            projector.load(target_aln_positions=target_aln_positions)

            for variant in block_variants:
                long_rows.extend(self.project_one_variant(variant, projector))

        return long_rows

    def project_one_variant(
        self,
        variant: VariantRecord,
        projector: AlignmentProjector,
    ) -> list[LongRow]:
        """Project one variant for all genotypes."""
        rows: list[LongRow] = []

        for genotype in self.genotype_order:
            projection = projector.get_projection(variant.aln_pos, genotype)
            block_start_in_region = self.block_starts_in_region.get((variant.block_id, genotype))
            region_start_in_source_seq = self.genotype_offsets.get(genotype, GenotypeOffset(genotype)).region_start_in_source_seq
            pos_in_region = compute_projected_position(block_start_in_region, projection.pos_in_block)
            pos_in_source_seq = compute_projected_position(
                region_start_in_source_seq,
                pos_in_region,
            )

            rows.append(
                LongRow(
                    block_id=variant.block_id,
                    aln_pos=variant.aln_pos,
                    genotype=genotype,
                    nt=projection.nt,
                    pos_in_block=projection.pos_in_block,
                    block_start_in_region=block_start_in_region,
                    pos_in_region=pos_in_region,
                    region_start_in_source_seq=region_start_in_source_seq,
                    pos_in_source_seq=pos_in_source_seq,
                )
            )

        return rows


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Project SNP positions from alignment coordinates to block, "
            "region, and source-sequence coordinates."
        )
    )
    parser.add_argument(
        "--vcf",
        required=True,
        type=Path,
        help="SNP VCF generated from block alignments.",
    )
    parser.add_argument(
        "--block-coords",
        required=True,
        type=Path,
        help="TSV with columns including: block_id, genotype, block_start_in_region.",
    )
    parser.add_argument(
        "--genotypes-tsv",
        required=True,
        type=Path,
        help=(
            "Input genotypes TSV used by the workflow. "
            "Column 2 must contain genotype names. Column 3 is optional and, if present, "
            "is interpreted as region_start_in_source_seq. Missing or empty offsets default to 1."
        ),
    )
    parser.add_argument(
        "--align-dir",
        required=True,
        type=Path,
        help="Directory containing per-block alignments named <block_id>.aln.fasta.",
    )
    parser.add_argument(
        "--long-output",
        required=True,
        type=Path,
        help="Output TSV in long format.",
    )
    parser.add_argument(
        "--wide-output",
        required=True,
        type=Path,
        help="Output TSV in wide format.",
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


def normalize_block_id(chrom_value: str) -> str:
    """Normalize a VCF CHROM value such as '4.aln' to the block identifier."""
    return chrom_value.removesuffix(".aln")


def read_vcf(vcf_path: Path) -> tuple[list[str], dict[str, list[VariantRecord]]]:
    """Read the VCF genotype order and group variants by block."""
    genotype_order: list[str] = []
    variants_by_block: dict[str, list[VariantRecord]] = {}

    with vcf_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("##"):
                continue

            if line.startswith("#CHROM"):
                fields = line.rstrip("\n").split("\t")
                genotype_order = fields[9:]
                continue

            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            chrom_value = fields[0]
            aln_pos = int(fields[1])
            block_id = normalize_block_id(chrom_value)
            variant = VariantRecord(block_id=block_id, aln_pos=aln_pos)
            variants_by_block.setdefault(block_id, []).append(variant)

    if not genotype_order:
        raise ValueError(f"Could not find VCF header with genotype names in {vcf_path}")

    LOGGER.info(
        "Read %d blocks and %d VCF genotypes from %s",
        len(variants_by_block),
        len(genotype_order),
        vcf_path,
    )
    return genotype_order, variants_by_block


def read_block_coords(block_coords_path: Path) -> dict[tuple[str, str], int]:
    """Read block start-in-region positions keyed by (block_id, genotype)."""
    dataframe = pl.read_csv(block_coords_path, separator="\t")
    required_columns = {"block_id", "genotype", "block_start_in_region"}

    if not required_columns.issubset(set(dataframe.columns)):
        raise ValueError(
            f"Missing required columns in block coordinates TSV {block_coords_path}: "
            f"expected {sorted(required_columns)}, got {dataframe.columns}"
        )

    block_starts_in_region: dict[tuple[str, str], int] = {}
    for row in dataframe.iter_rows(named=True):
        key = (str(row["block_id"]), str(row["genotype"]))
        block_starts_in_region[key] = int(row["block_start_in_region"])

    LOGGER.info(
        "Read %d block coordinate entries from %s",
        len(block_starts_in_region),
        block_coords_path,
    )
    return block_starts_in_region


def read_genotype_offsets(genotypes_tsv_path: Path) -> dict[str, GenotypeOffset]:
    """Read per-genotype source-sequence offsets from the workflow genotypes TSV."""
    genotype_offsets: dict[str, GenotypeOffset] = {}

    with genotypes_tsv_path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            genotype = (row.get("genotype") or "").strip()
            region_start_in_source_seq = int(row["region_start"]) if (row.get("region_start") or "").strip() else 1

            genotype_offsets[genotype] = GenotypeOffset(
                genotype=genotype,
                region_start_in_source_seq=region_start_in_source_seq,
            )

    LOGGER.info("Read %d genotype offsets from %s", len(genotype_offsets), genotypes_tsv_path)
    return genotype_offsets


def compute_projected_position(start_position: int | None, relative_position: int | None) -> int | None:
    """Project a relative position onto a coordinate system using a 1-based start."""
    if start_position is None or relative_position is None:
        return None
    return start_position + relative_position - 1


def natural_sort_key(value: str) -> tuple[int, str]:
    """Sort numerically when possible, otherwise lexicographically."""
    if value.isdigit():
        return int(value), value
    return 10**18, value


def build_long_dataframe(long_rows: list[LongRow]) -> pl.DataFrame:
    """Build the long-format output dataframe."""
    rows = [
        {
            "snp_id": row.snp_id,
            "block_id": row.block_id,
            "aln_pos": row.aln_pos,
            "genotype": row.genotype,
            "nt": row.nt,
            "pos_in_block": row.pos_in_block,
            "block_start_in_region": row.block_start_in_region,
            "pos_in_region": row.pos_in_region,
            "region_start_in_source_seq": row.region_start_in_source_seq,
            "pos_in_source_seq": row.pos_in_source_seq,
        }
        for row in long_rows
    ]

    dataframe = pl.DataFrame(rows)
    return dataframe.sort(["block_id", "aln_pos", "genotype"])


def build_wide_dataframe(long_rows: list[LongRow], genotype_order: list[str]) -> pl.DataFrame:
    """Build the wide-format output dataframe."""
    grouped_rows: dict[tuple[str, int], dict[str, str | int | None]] = {}

    for row in long_rows:
        key = (row.block_id, row.aln_pos)

        if key not in grouped_rows:
            grouped_rows[key] = {
                "snp_id": row.snp_id,
                "block_id": row.block_id,
                "aln_pos": row.aln_pos,
            }

        grouped_rows[key][f"{row.genotype}_nt"] = row.nt
        grouped_rows[key][f"{row.genotype}_pos"] = row.pos_in_source_seq

    wide_rows = [
        grouped_rows[key]
        for key in sorted(grouped_rows, key=lambda x: (natural_sort_key(x[0]), x[1]))
    ]
    dataframe = pl.DataFrame(wide_rows)

    ordered_columns = ["snp_id", "block_id", "aln_pos"]
    ordered_columns.extend(f"{genotype}_nt" for genotype in genotype_order)
    ordered_columns.extend(f"{genotype}_pos" for genotype in genotype_order)

    return dataframe.select(ordered_columns)


def write_dataframe(dataframe: pl.DataFrame, output_path: Path) -> None:
    """Write a dataframe as a TSV file."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    dataframe.write_csv(output_path, separator="\t")
    LOGGER.info("Wrote %s", output_path)


def main() -> None:
    """Run the SNP position projection script."""
    args = parse_args()
    setup_logging(args.log_level)

    projector = SnpPositionProjector(
        vcf_path=args.vcf,
        block_coords_path=args.block_coords,
        genotypes_tsv_path=args.genotypes_tsv,
        align_dir=args.align_dir,
        long_output_path=args.long_output,
        wide_output_path=args.wide_output,
    )
    projector.run()


if __name__ == "__main__":
    main()
