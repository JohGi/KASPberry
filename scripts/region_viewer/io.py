#!/usr/bin/env python3
# Author: Johanna Girodolle

from __future__ import annotations

from pathlib import Path

import polars as pl
import json
import logging
import re
import os
import csv
from urllib.parse import unquote

from .models import BlockAlignment, BlockFeature, DistanceMatrix, DotplotRecord, GffGeneFeature, GffTrack, SampleData, SampleRecord, SnpFeature

LOGGER = logging.getLogger(__name__)

COORD_SUFFIX_PATTERN = re.compile(r":\d+-\d+$")


def read_samples(path: Path) -> list[SampleRecord]:
    """Read genotype definitions from the KASPberry genotype table."""
    records: list[SampleRecord] = []

    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        if reader.fieldnames is None:
            raise ValueError(f"Genotype table has no header: {path}")

        required_columns = {"genotype", "region_fasta", "region_start"}
        missing_columns = required_columns - set(reader.fieldnames)

        if missing_columns:
            raise ValueError(
                f"Missing required columns in genotype table {path}: "
                f"{sorted(missing_columns)}"
            )

        for row in reader:
            records.append(
                SampleRecord(
                    fasta_path=Path(row["region_fasta"]),
                    sample=row["genotype"],
                    region_start_in_source_seq=(
                        int(row["region_start"])
                        if (row.get("region_start") or "").strip()
                        else 1
                    ),
                )
            )

    if not records:
        raise ValueError(f"Genotype table is empty: {path}")

    return records


def read_single_fasta_length(path: Path) -> int:
    """Read the sequence length of a single-sequence FASTA."""
    if not path.is_file():
        raise FileNotFoundError(f"FASTA file not found: {path}")

    header_count = 0
    sequence_chunks: list[str] = []

    with path.open(encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith(">"):
                header_count += 1
                continue
            sequence_chunks.append(stripped)

    if header_count != 1:
        raise ValueError(f"Expected exactly one sequence in FASTA {path}, found {header_count}")

    sequence = "".join(sequence_chunks)
    if not sequence:
        raise ValueError(f"Empty sequence in FASTA: {path}")

    return len(sequence)


def read_fasta_lengths(fasta_dir: Path) -> dict[str, int]:
    """Read lengths from per-sample FASTA files."""
    if not fasta_dir.is_dir():
        raise FileNotFoundError(f"FASTA directory not found: {fasta_dir}")

    lengths: dict[str, int] = {}
    for fasta_path in sorted(fasta_dir.glob("*.fasta")):
        lengths[fasta_path.stem] = read_single_fasta_length(fasta_path)

    if not lengths:
        raise ValueError(f"No FASTA files found in directory: {fasta_dir}")

    return lengths


def read_blocks(path: Path) -> dict[str, list[BlockFeature]]:
    """Read collinear blocks from the enriched block coordinates TSV."""
    dataframe = pl.read_csv(path, separator="\t")

    required_columns = {
        "block_id",
        "genotype",
        "block_start_in_region",
        "block_end_in_region",
        "block_start_in_source_seq",
        "block_end_in_source_seq",
    }
    missing_columns = required_columns - set(dataframe.columns)
    if missing_columns:
        raise ValueError(
            f"Missing required columns in block coordinates TSV {path}: {sorted(missing_columns)}"
        )

    blocks_by_sample: dict[str, list[BlockFeature]] = {}
    for row in dataframe.iter_rows(named=True):
        block = BlockFeature(
            sample=str(row["genotype"]),
            block_id=str(row["block_id"]),
            block_start_in_region=int(row["block_start_in_region"]),
            block_end_in_region=int(row["block_end_in_region"]),
            block_start_in_source_seq=int(row["block_start_in_source_seq"]),
            block_end_in_source_seq=int(row["block_end_in_source_seq"]),
        )
        blocks_by_sample.setdefault(block.sample, []).append(block)

    return blocks_by_sample


def read_snp_long(path: Path) -> pl.DataFrame:
    """Load and validate the long-format SNP TSV as a polars DataFrame."""
    dataframe = pl.read_csv(path, separator="\t")
    required_columns = {
        "block_id",
        "aln_pos",
        "genotype",
        "nt",
        "pos_in_block",
        "pos_in_region",
        "pos_in_source_seq",
    }
    missing_columns = required_columns - set(dataframe.columns)
    if missing_columns:
        raise ValueError(f"Missing required columns in SNP TSV {path}: {sorted(missing_columns)}")
    return dataframe


def parse_snps(snp_long: pl.DataFrame) -> dict[str, list[SnpFeature]]:
    """Convert a long-format SNP DataFrame into per-sample SnpFeature lists."""
    snps_by_sample: dict[str, list[SnpFeature]] = {}
    for row in snp_long.iter_rows(named=True):
        snp = SnpFeature(
            sample=str(row["genotype"]),
            block_id=str(row["block_id"]),
            aln_pos=int(row["aln_pos"]),
            nt=str(row["nt"]),
            pos_in_block=int(row["pos_in_block"]),
            pos_in_region=int(row["pos_in_region"]),
            pos_in_source_seq=int(row["pos_in_source_seq"]),
        )
        snps_by_sample.setdefault(snp.sample, []).append(snp)
    return snps_by_sample


def count_unique_snps(snp_long: pl.DataFrame) -> int:
    """Count unique SNP markers identified by (block_id, aln_pos) pairs."""
    if snp_long.is_empty():
        return 0
    return snp_long.select(["block_id", "aln_pos"]).unique().height


def read_snps(path: Path) -> dict[str, list[SnpFeature]]:
    """Read SNPs from the long-format SNP TSV."""
    return parse_snps(read_snp_long(path))


def read_summary_stats(path: Path) -> dict[str, object]:
    """Read summary statistics from a JSON file."""
    return json.loads(path.read_text(encoding="utf-8"))


def read_masked_block_n_stats(
    path: Path,
) -> dict[str, dict[str, dict[str, int | float]]]:
    """Read masked block N-content statistics as a nested dictionary."""
    dataframe = pl.read_csv(path, separator="\t")

    stats: dict[str, dict[str, dict[str, int | float]]] = {}

    for row in dataframe.iter_rows(named=True):
        block_id = str(row["block_id"])
        sample = str(row["sample"])

        stats.setdefault(block_id, {})[sample] = {
            "length_bp": int(row["masked_length_bp"]),
            "unmasked_n_count": int(row["unmasked_n_count"]),
            "unmasked_n_pct": float(row["unmasked_n_pct"]),
            "masked_n_count": int(row["masked_n_count"]),
            "masked_n_pct": float(row["masked_n_pct"]),
            "repeat_masked_n_count": int(row["repeat_masked_n_count"]),
            "repeat_masked_n_pct": float(row["repeat_masked_n_pct"]),
        }

    return stats


def normalize_distance(value: float) -> float:
    """Normalize near-zero distances."""
    return 0.0 if abs(value) < 1e-9 else value


def parse_mash_matrix(path: Path, sample_order: list[str]) -> DistanceMatrix:
    """Parse a Mash square distance matrix."""
    rows = path.read_text(encoding="utf-8").splitlines()

    header = [
        Path(label).stem
        for label in rows[0].split("\t")[1:]
    ]

    values_by_pair: dict[tuple[str, str], float] = {}

    for row in rows[1:]:
        fields = row.split("\t")
        row_label = Path(fields[0]).stem
        distances = [float(value) for value in fields[1:]]

        for col_label, distance in zip(header, distances):
            values_by_pair[(row_label, col_label)] = normalize_distance(distance)

    values = [
        [
            values_by_pair[(row_label, col_label)]
            for col_label in sample_order
        ]
        for row_label in sample_order
    ]

    return DistanceMatrix(
        labels=sample_order,
        values=values,
        source="mash",
        title="Mash distances, whole region",
        unit="mash_distance",
    )


def parse_emboss_distmat(path: Path, sample_order: list[str], block_id: str) -> DistanceMatrix:
    """Parse an EMBOSS distmat triangular matrix."""
    lines = path.read_text(encoding="utf-8").splitlines()

    labels: list[str] = []
    triangle_rows: list[list[float]] = []

    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue

        tokens = stripped.split()
        if len(tokens) < 2:
            continue

        if not tokens[-1].isdigit():
            continue

        label = tokens[-2]
        numeric_part = tokens[:-2]

        try:
            distances = [
                normalize_distance(float(value) / 100)
                for value in numeric_part
            ]
        except ValueError:
            continue

        labels.append(label)
        triangle_rows.append(distances)

    n = len(labels)
    matrix_by_label = {
        label: {other_label: 0.0 for other_label in labels}
        for label in labels
    }

    for row_index, row_values in enumerate(triangle_rows):
        row_label = labels[row_index]

        for offset, distance in enumerate(row_values):
            col_index = row_index + offset
            col_label = labels[col_index]

            matrix_by_label[row_label][col_label] = distance
            matrix_by_label[col_label][row_label] = distance

    values = [
        [
            matrix_by_label[row_label][col_label]
            for col_label in sample_order
        ]
        for row_label in sample_order
    ]

    return DistanceMatrix(
        labels=sample_order,
        values=values,
        source="kimura2p",
        title="Kimura 2P distances",
        unit="substitutions_per_base",
    )


def parse_kimura2p_distmat_dir(
    distmat_dir: Path,
    sample_order: list[str],
) -> dict[str, dict[str, object]]:
    """Parse all Kimura 2P EMBOSS distmat files from a directory."""
    matrices = {}

    for path in sorted(distmat_dir.glob("*.kimura2p.distmat")):
        block_id = path.name.split(".", 1)[0]
        matrices[block_id] = parse_emboss_distmat(
            path=path,
            sample_order=sample_order,
            block_id=block_id,
        ).to_dict()

    return matrices

def parse_alignment_sample_name(header: str) -> str:
    """Extract the sample name from an alignment FASTA header."""
    first_token = header.strip().split()[0]
    return COORD_SUFFIX_PATTERN.sub("", first_token)


def read_fasta_alignment(path: Path) -> dict[str, str]:
    """Read a FASTA alignment as a sample-to-sequence mapping."""
    if not path.is_file():
        raise FileNotFoundError(f"Alignment file not found: {path}")

    sequences: dict[str, str] = {}
    current_name: str | None = None
    current_chunks: list[str] = []

    with path.open(encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()

            if not stripped:
                continue

            if stripped.startswith(">"):
                if current_name is not None:
                    sequences[current_name] = "".join(current_chunks)

                current_name = parse_alignment_sample_name(stripped[1:])
                current_chunks = []
                continue

            current_chunks.append(stripped)

    if current_name is not None:
        sequences[current_name] = "".join(current_chunks)

    if not sequences:
        raise ValueError(f"No aligned sequences found in alignment file: {path}")

    return sequences


def read_block_alignments(
    align_dir: Path,
    block_ids: list[str],
) -> dict[str, BlockAlignment]:
    """Read masked block alignments indexed by block ID."""
    if not align_dir.is_dir():
        raise FileNotFoundError(f"Alignment directory not found: {align_dir}")

    alignments: dict[str, BlockAlignment] = {}

    for block_id in block_ids:
        alignment_path = align_dir / f"{block_id}.aln.fasta"

        if not alignment_path.exists():
            LOGGER.warning("Missing alignment for block %s: %s", block_id, alignment_path)
            continue

        alignments[block_id] = BlockAlignment(
            block_id=block_id,
            sequences_by_sample=read_fasta_alignment(alignment_path),
        )

    return alignments


def read_gff_tracks_json(path: Path) -> dict[str, dict[str, Path]]:
    """Read configured GFF tracks from a JSON file."""
    if not path.is_file():
        raise FileNotFoundError(f"GFF tracks JSON file not found: {path}")

    with path.open(encoding="utf-8") as handle:
        raw_tracks = json.load(handle)

    if raw_tracks is None:
        return {}

    if not isinstance(raw_tracks, dict):
        raise ValueError("GFF tracks JSON content must be a dictionary.")

    gff_tracks: dict[str, dict[str, Path]] = {}

    for sample_name, sample_tracks in raw_tracks.items():
        if not isinstance(sample_name, str):
            raise ValueError("GFF track sample names must be strings.")

        if not isinstance(sample_tracks, dict):
            raise ValueError(
                f"GFF tracks for sample {sample_name!r} must be a dictionary."
            )

        gff_tracks[sample_name] = {}

        for track_name, gff_path in sample_tracks.items():
            if not isinstance(track_name, str):
                raise ValueError(
                    f"GFF track names for sample {sample_name!r} must be strings."
                )

            if not isinstance(gff_path, str):
                raise ValueError(
                    f"GFF path for {sample_name}.{track_name} must be a string."
                )

            gff_tracks[sample_name][track_name] = Path(gff_path)

    return gff_tracks


def parse_gff_attributes(attributes: str) -> dict[str, str]:
    """Parse a GFF3 attribute field into a dictionary."""
    parsed_attributes: dict[str, str] = {}

    if attributes == ".":
        return parsed_attributes

    for item in attributes.split(";"):
        if not item:
            continue

        if "=" not in item:
            parsed_attributes[unquote(item)] = ""
            continue

        key, value = item.split("=", 1)
        parsed_attributes[unquote(key)] = unquote(value)

    return parsed_attributes


def get_gene_id(attributes: dict[str, str], fallback_id: str) -> str:
    """Return the best available gene identifier from GFF attributes."""
    for key in ("ID", "Name", "gene_id", "locus_tag"):
        value = attributes.get(key)

        if value:
            return value

    return fallback_id


def overlaps_region(
    feature_start: int,
    feature_end: int,
    region_start: int,
    region_end: int,
) -> bool:
    """Return whether a source-sequence feature overlaps a source-sequence region."""
    return feature_end >= region_start and feature_start <= region_end


def project_source_interval_to_region(
    source_start: int,
    source_end: int,
    region_start: int,
    region_length: int,
) -> tuple[int, int]:
    """Project and clip a source-sequence interval into region coordinates."""
    raw_start_in_region = source_start - region_start + 1
    raw_end_in_region = source_end - region_start + 1

    clipped_start_in_region = max(1, raw_start_in_region)
    clipped_end_in_region = min(region_length, raw_end_in_region)

    return clipped_start_in_region, clipped_end_in_region


def read_projected_gff_gene_features(
    path: Path,
    sample: str,
    track_name: str,
    region_start_in_source_seq: int,
    region_length: int,
) -> GffTrack:
    """Read one GFF file and project gene features into region coordinates."""
    if not path.is_file():
        raise FileNotFoundError(f"GFF file not found: {path}")

    region_end_in_source_seq = region_start_in_source_seq + region_length - 1
    source_seq_ids: set[str] = set()
    features: list[GffGeneFeature] = []

    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.strip()

            if not stripped or stripped.startswith("#"):
                continue

            fields = stripped.split("\t")
            if len(fields) != 9:
                raise ValueError(
                    f"Invalid GFF line in {path} at line {line_number}: "
                    f"expected 9 columns, got {len(fields)}."
                )

            (
                source_seq_id,
                _source,
                feature_type,
                start,
                end,
                _score,
                strand,
                _phase,
                attributes,
            ) = fields

            source_seq_ids.add(source_seq_id)

            if feature_type != "gene":
                continue

            start_in_source_seq = int(start)
            end_in_source_seq = int(end)

            if start_in_source_seq > end_in_source_seq:
                raise ValueError(
                    f"Invalid gene coordinates in {path} at line {line_number}: "
                    f"start {start_in_source_seq} is greater than end {end_in_source_seq}."
                )

            if not overlaps_region(
                feature_start=start_in_source_seq,
                feature_end=end_in_source_seq,
                region_start=region_start_in_source_seq,
                region_end=region_end_in_source_seq,
            ):
                continue

            start_in_region, end_in_region = project_source_interval_to_region(
                source_start=start_in_source_seq,
                source_end=end_in_source_seq,
                region_start=region_start_in_source_seq,
                region_length=region_length,
            )

            parsed_attributes = parse_gff_attributes(attributes)
            fallback_id = f"{sample}_{track_name}_gene_{line_number}"

            features.append(
                GffGeneFeature(
                    sample=sample,
                    track_name=track_name,
                    gene_id=get_gene_id(parsed_attributes, fallback_id),
                    source_seq_id=source_seq_id,
                    start_in_source_seq=start_in_source_seq,
                    end_in_source_seq=end_in_source_seq,
                    start_in_region=start_in_region,
                    end_in_region=end_in_region,
                    strand=None if strand == "." else strand,
                )
            )

    if len(source_seq_ids) > 1:
        raise ValueError(
            f"GFF file must contain a single sequence ID, but {path} contains: "
            + ", ".join(sorted(source_seq_ids))
        )

    return GffTrack(
        sample=sample,
        track_name=track_name,
        features=features,
    )


def read_gff_gene_tracks(
    gff_tracks: dict[str, dict[str, Path]],
    sample_data: list[SampleData],
) -> dict[str, list[GffTrack]]:
    """Read and project configured GFF gene tracks by sample."""
    samples_by_name = {sample.sample: sample for sample in sample_data}
    projected_tracks: dict[str, list[GffTrack]] = {
        sample.sample: []
        for sample in sample_data
    }

    for sample_name, sample_tracks in gff_tracks.items():
        if sample_name not in samples_by_name:
            raise ValueError(f"Unknown sample in GFF tracks JSON: {sample_name}")

        sample = samples_by_name[sample_name]

        for track_name, gff_path in sample_tracks.items():
            projected_tracks[sample_name].append(
                read_projected_gff_gene_features(
                    path=gff_path,
                    sample=sample_name,
                    track_name=track_name,
                    region_start_in_source_seq=sample.region_start_in_source_seq,
                    region_length=sample.region_length,
                )
            )

    return projected_tracks


def read_dotplot_manifest(path: Path | None, output_path: Path) -> list[DotplotRecord]:
    """Read dotplot records from a manifest JSON, resolving SVG paths relative to output_path."""
    if path is None:
        return []

    if not path.is_file():
        raise FileNotFoundError(f"Dotplot manifest JSON not found: {path}")

    raw = json.loads(path.read_text(encoding="utf-8"))

    if not isinstance(raw, dict) or "dotplots" not in raw:
        raise ValueError("Dotplot manifest JSON must be a root object with a 'dotplots' list.")

    raw_pairs = raw["dotplots"]
    if not isinstance(raw_pairs, list):
        raise ValueError("The 'dotplots' key in the manifest JSON must be a list.")

    required_keys = {"pair_id", "x_sample", "y_sample", "svg_rel_path"}
    records: list[DotplotRecord] = []

    for index, entry in enumerate(raw_pairs):
        missing = required_keys - set(entry.keys())
        if missing:
            raise ValueError(
                f"Dotplot entry at index {index} is missing required keys: {sorted(missing)}"
            )

        absolute_svg = (path.parent / entry["svg_rel_path"]).resolve()
        svg_rel_path = os.path.relpath(
            absolute_svg,
            start=output_path.parent.resolve(),
        )

        records.append(
            DotplotRecord(
                pair_id=str(entry["pair_id"]),
                x_sample=str(entry["x_sample"]),
                y_sample=str(entry["y_sample"]),
                svg_rel_path=svg_rel_path,
            )
        )

    return records


def write_html(html: str, output_path: Path) -> None:
    """Write the final HTML file."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(html, encoding="utf-8")
    LOGGER.info("Wrote %s", output_path)
