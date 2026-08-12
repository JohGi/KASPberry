#!/usr/bin/env python3
# Author: Johanna Girodolle

"""Map filtered SNP coordinates from alignment coordinates to block, zone, and source sequence coordinates."""

from __future__ import annotations

import argparse
import csv
import logging
from pathlib import Path

import polars as pl
from attrs import define, field
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

LOGGER = logging.getLogger(__name__)


@define(frozen=True)
class SampleOffset:
    """Store the source-sequence start offset for one sample."""

    sample: str
    zone_start_in_source_seq: int = 1


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
    sample: str
    nt: str
    pos_in_block: int | None
    block_start_in_zone: int | None
    pos_in_zone: int | None
    zone_start_in_source_seq: int | None
    pos_in_source_seq: int | None


@define(frozen=True)
class SnpProjection:
    """Store projected values for one sample at one SNP position."""

    nt: str
    pos_in_block: int | None


@define
class AlignmentProjector:
    """Project selected alignment columns to ungapped positions for each sample."""

    block_id: str
    alignment_path: Path
    sample_order: list[str]
    projections_by_aln_pos: dict[int, dict[str, SnpProjection]] = field(factory=dict)

    def load(self, target_aln_positions: set[int]) -> None:
        """Load the alignment and precompute projections only for requested SNP columns."""
        if not target_aln_positions:
            self.projections_by_aln_pos = {}
            return

        alignment = AlignIO.read(str(self.alignment_path), "fasta")
        normalized_names = self.get_normalized_sample_names(alignment)
        self.validate_alignment_samples(normalized_names)
        self.projections_by_aln_pos = self.build_projection_cache(
            alignment=alignment,
            normalized_names=normalized_names,
            target_aln_positions=target_aln_positions,
        )

    @staticmethod
    def normalize_alignment_sample_name(sequence_id: str) -> str:
        """Extract the sample name from an alignment record identifier."""
        return sequence_id.split(":", 1)[0]

    def get_normalized_sample_names(self, alignment: MultipleSeqAlignment) -> list[str]:
        """Return normalized sample names for all alignment records."""
        return [self.normalize_alignment_sample_name(record.id) for record in alignment]

    def validate_alignment_samples(self, normalized_names: list[str]) -> None:
        """Validate alignment sample names against the expected VCF sample order."""
        observed_samples = set(normalized_names)
        expected_samples = set(self.sample_order)

        if len(observed_samples) != len(normalized_names):
            raise ValueError(
                f"Duplicate normalized sample names found in alignment for block {self.block_id}: "
                f"{normalized_names}"
            )

        if observed_samples != expected_samples:
            raise ValueError(
                f"Alignment samples do not match VCF samples for block {self.block_id}. "
                f"Expected={sorted(expected_samples)} Observed={sorted(observed_samples)}"
            )

    def build_projection_cache(
        self,
        alignment: MultipleSeqAlignment,
        normalized_names: list[str],
        target_aln_positions: set[int],
    ) -> dict[int, dict[str, SnpProjection]]:
        """Build per-sample projections only for requested alignment positions."""
        alignment_length = alignment.get_alignment_length()
        max_target = max(target_aln_positions)
        ungapped_counters: dict[str, int] = {sample: 0 for sample in normalized_names}
        projections_by_aln_pos: dict[int, dict[str, SnpProjection]] = {}

        for aln_index in range(alignment_length):
            aln_pos = aln_index + 1

            for record, sample in zip(alignment, normalized_names):
                nt = str(record.seq[aln_index]).upper()
                if nt != "-":
                    ungapped_counters[sample] += 1

            if aln_pos not in target_aln_positions:
                if aln_pos >= max_target:
                    break
                continue

            projections_by_aln_pos[aln_pos] = {}

            for record, sample in zip(alignment, normalized_names):
                nt = str(record.seq[aln_index]).upper()
                pos_in_block = ungapped_counters[sample] if nt != "-" else None
                projections_by_aln_pos[aln_pos][sample] = SnpProjection(
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

    def get_projection(self, aln_pos: int, sample: str) -> SnpProjection:
        """Return the nucleotide and ungapped block position for one sample at one SNP column."""
        return self.projections_by_aln_pos[aln_pos][sample]


@define
class SnpPositionProjector:
    """Project SNP coordinates from alignments to zone and source sequences."""

    vcf_path: Path
    block_coords_path: Path
    samples_tsv_path: Path
    align_dir: Path
    long_output_path: Path
    wide_output_path: Path
    sample_order: list[str] = field(factory=list)
    variants_by_block: dict[str, list[VariantRecord]] = field(factory=dict)
    block_starts_in_zone: dict[tuple[str, str], int] = field(factory=dict)
    sample_offsets: dict[str, SampleOffset] = field(factory=dict)

    def run(self) -> None:
        """Run the full SNP projection workflow."""
        self.sample_order, self.variants_by_block = read_vcf(self.vcf_path)
        self.block_starts_in_zone = read_block_coords(self.block_coords_path)
        self.sample_offsets = read_sample_offsets(self.samples_tsv_path)
        long_rows = self.project_variants()
        long_df = build_long_dataframe(long_rows)
        wide_df = build_wide_dataframe(long_rows, self.sample_order)
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
                sample_order=self.sample_order,
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
        """Project one variant for all samples."""
        rows: list[LongRow] = []

        for sample in self.sample_order:
            projection = projector.get_projection(variant.aln_pos, sample)
            block_start_in_zone = self.block_starts_in_zone.get((variant.block_id, sample))
            zone_start_in_source_seq = self.sample_offsets.get(sample, SampleOffset(sample)).zone_start_in_source_seq
            pos_in_zone = compute_projected_position(block_start_in_zone, projection.pos_in_block)
            pos_in_source_seq = compute_projected_position(
                zone_start_in_source_seq,
                pos_in_zone,
            )

            rows.append(
                LongRow(
                    block_id=variant.block_id,
                    aln_pos=variant.aln_pos,
                    sample=sample,
                    nt=projection.nt,
                    pos_in_block=projection.pos_in_block,
                    block_start_in_zone=block_start_in_zone,
                    pos_in_zone=pos_in_zone,
                    zone_start_in_source_seq=zone_start_in_source_seq,
                    pos_in_source_seq=pos_in_source_seq,
                )
            )

        return rows


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Project filtered SNP positions from alignment coordinates to block, "
            "zone, and source-sequence coordinates."
        )
    )
    parser.add_argument(
        "--vcf",
        required=True,
        type=Path,
        help="Filtered or unfiltered SNP VCF generated from block alignments.",
    )
    parser.add_argument(
        "--block-coords",
        required=True,
        type=Path,
        help="TSV with columns including: block_id, sample, block_start_in_zone.",
    )
    parser.add_argument(
        "--samples-tsv",
        required=True,
        type=Path,
        help=(
            "Input samples TSV used by the workflow. "
            "Column 2 must contain sample names. Column 3 is optional and, if present, "
            "is interpreted as zone_start_in_source_seq. Missing or empty offsets default to 1."
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
    """Read the VCF sample order and group variants by block."""
    sample_order: list[str] = []
    variants_by_block: dict[str, list[VariantRecord]] = {}

    with vcf_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("##"):
                continue

            if line.startswith("#CHROM"):
                fields = line.rstrip("\n").split("\t")
                sample_order = fields[9:]
                continue

            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            chrom_value = fields[0]
            aln_pos = int(fields[1])
            block_id = normalize_block_id(chrom_value)
            variant = VariantRecord(block_id=block_id, aln_pos=aln_pos)
            variants_by_block.setdefault(block_id, []).append(variant)

    if not sample_order:
        raise ValueError(f"Could not find VCF header with sample names in {vcf_path}")

    LOGGER.info(
        "Read %d blocks and %d VCF samples from %s",
        len(variants_by_block),
        len(sample_order),
        vcf_path,
    )
    return sample_order, variants_by_block


def read_block_coords(block_coords_path: Path) -> dict[tuple[str, str], int]:
    """Read block start-in-zone positions keyed by (block_id, sample)."""
    dataframe = pl.read_csv(block_coords_path, separator="\t")
    required_columns = {"block_id", "sample", "block_start_in_zone"}

    if not required_columns.issubset(set(dataframe.columns)):
        raise ValueError(
            f"Missing required columns in block coordinates TSV {block_coords_path}: "
            f"expected {sorted(required_columns)}, got {dataframe.columns}"
        )

    block_starts_in_zone: dict[tuple[str, str], int] = {}
    for row in dataframe.iter_rows(named=True):
        key = (str(row["block_id"]), str(row["sample"]))
        block_starts_in_zone[key] = int(row["block_start_in_zone"])

    LOGGER.info(
        "Read %d block coordinate entries from %s",
        len(block_starts_in_zone),
        block_coords_path,
    )
    return block_starts_in_zone


def read_sample_offsets(samples_tsv_path: Path) -> dict[str, SampleOffset]:
    """Read per-sample source-sequence offsets from the workflow samples TSV."""
    sample_offsets: dict[str, SampleOffset] = {}

    with samples_tsv_path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            sample = (row.get("genotype") or "").strip()
            zone_start_in_source_seq = int(row["region_start"]) if (row.get("region_start") or "").strip() else 1

            sample_offsets[sample] = SampleOffset(
                sample=sample,
                zone_start_in_source_seq=zone_start_in_source_seq,
            )

    LOGGER.info("Read %d sample offsets from %s", len(sample_offsets), samples_tsv_path)
    return sample_offsets


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
            "block_id": row.block_id,
            "aln_pos": row.aln_pos,
            "sample": row.sample,
            "nt": row.nt,
            "pos_in_block": row.pos_in_block,
            "block_start_in_zone": row.block_start_in_zone,
            "pos_in_zone": row.pos_in_zone,
            "zone_start_in_source_seq": row.zone_start_in_source_seq,
            "pos_in_source_seq": row.pos_in_source_seq,
        }
        for row in long_rows
    ]

    dataframe = pl.DataFrame(rows)
    return dataframe.sort(["block_id", "aln_pos", "sample"])


def build_wide_dataframe(long_rows: list[LongRow], sample_order: list[str]) -> pl.DataFrame:
    """Build the wide-format output dataframe."""
    grouped_rows: dict[tuple[str, int], dict[str, str | int | None]] = {}

    for row in long_rows:
        key = (row.block_id, row.aln_pos)

        if key not in grouped_rows:
            grouped_rows[key] = {
                "block_id": row.block_id,
                "aln_pos": row.aln_pos,
            }

        grouped_rows[key][f"{row.sample}_nt"] = row.nt
        grouped_rows[key][f"{row.sample}_pos"] = row.pos_in_source_seq

    wide_rows = [
        grouped_rows[key]
        for key in sorted(grouped_rows, key=lambda x: (natural_sort_key(x[0]), x[1]))
    ]
    dataframe = pl.DataFrame(wide_rows)

    ordered_columns = ["block_id", "aln_pos"]
    ordered_columns.extend(f"{sample}_nt" for sample in sample_order)
    ordered_columns.extend(f"{sample}_pos" for sample in sample_order)

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
        samples_tsv_path=args.samples_tsv,
        align_dir=args.align_dir,
        long_output_path=args.long_output,
        wide_output_path=args.wide_output,
    )
    projector.run()


if __name__ == "__main__":
    main()
