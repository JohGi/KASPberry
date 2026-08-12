#!/usr/bin/env python3
# Author: Johanna Girodolle

"""Build a JSON manifest for pairwise dotplot-only SVG files."""

from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Build a JSON manifest for pairwise dotplot-only SVG files."
    )
    parser.add_argument("--genotypes", required=True, help="Input genotype TSV.")
    parser.add_argument("--svg-dir", required=True, help="Directory containing SVG files.")
    parser.add_argument("--output", required=True, help="Output JSON manifest path.")
    return parser.parse_args()


def read_genotype_names(genotypes_path: Path) -> list[str]:
    """Read genotype names from the KASPberry genotype table."""
    genotype_names: list[str] = []

    with genotypes_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        if reader.fieldnames is None:
            raise ValueError(f"Genotype table has no header: {genotypes_path}")

        if "genotype" not in reader.fieldnames:
            raise ValueError(
                f"Missing required 'genotype' column in {genotypes_path}"
            )

        for row in reader:
            genotype = row["genotype"].strip()

            if not genotype:
                raise ValueError(
                    f"Empty genotype name in {genotypes_path}"
                )

            genotype_names.append(genotype)

    if not genotype_names:
        raise ValueError(f"Genotype table is empty: {genotypes_path}")

    return genotype_names


def build_pair_id(sample_a: str, sample_b: str) -> str:
    """Build the pair identifier used by the dotplot workflow."""
    return f"{sample_a}__vs__{sample_b}"


def build_dotplot_records(
    sample_names: list[str],
    svg_dir: Path,
    output_path: Path,
) -> list[dict[str, str]]:
    """Build manifest records for existing dotplot-only SVG files."""
    records: list[dict[str, str]] = []

    for row_index, row_sample in enumerate(sample_names[:-1]):
        for col_sample in sample_names[row_index + 1:]:
            pair_id = build_pair_id(row_sample, col_sample)
            svg_path = svg_dir / f"{pair_id}.dotplot_only.svg"

            if not svg_path.exists():
                continue

            records.append(
                {
                    "pair_id": pair_id,
                    "x_sample": col_sample,
                    "y_sample": row_sample,
                    "svg_rel_path": os.path.relpath(svg_path, start=output_path.parent),
                }
            )

    return records


def write_manifest(output_path: Path, records: list[dict[str, str]]) -> None:
    """Write the dotplot manifest as JSON."""
    payload = {
        "format_version": 1,
        "dotplots": records,
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )


def main() -> None:
    """Run the manifest builder."""
    args = parse_args()
    output_path = Path(args.output)

    genotype_names = read_genotype_names(Path(args.genotypes))

    records = build_dotplot_records(
        sample_names=genotype_names,
        svg_dir=Path(args.svg_dir),
        output_path=output_path,
    )
    write_manifest(output_path, records)


if __name__ == "__main__":
    main()
