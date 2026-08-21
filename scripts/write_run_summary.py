#!/usr/bin/env python3
# Author: Johanna Girodolle

"""Write a concise summary of a KASPberry run."""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

import polars as pl


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--mode",
        choices=["snps", "kasp"],
        required=True,
    )

    parser.add_argument(
        "--repeat-masking",
        action="store_true",
    )
    parser.add_argument(
        "--repeat-masking-library",
        type=Path,
    )

    parser.add_argument(
        "--genotypes",
        type=Path,
        required=True,
    )
    parser.add_argument(
        "--block-coords",
        type=Path,
        required=True,
    )
    parser.add_argument(
        "--snp-summary",
        type=Path,
        required=True,
    )
    parser.add_argument(
        "--assay-summary",
        type=Path,
    )
    parser.add_argument(
        "--masked-block-n-stats",
        type=Path,
        required=True,
    )
    parser.add_argument(
        "--clean-fastas",
        type=Path,
        nargs="+",
        required=True,
    )
    parser.add_argument(
        "--json-output",
        type=Path,
        required=True,
    )
    parser.add_argument(
        "--txt-output",
        type=Path,
        required=True,
    )

    return parser.parse_args()


def fasta_length(path: Path) -> int:
    return sum(
        len(line.strip())
        for line in path.read_text().splitlines()
        if line and not line.startswith(">")
    )


def count_reasons(
    df: pl.DataFrame,
    status_col: str,
    reason_col: str,
) -> dict[str, int]:
    reasons = Counter(
        row[reason_col]
        for row in df.iter_rows(named=True)
        if row[status_col] == "FAIL" and row[reason_col]
    )

    return dict(sorted(reasons.items()))


def interval_union_length(
    intervals: list[tuple[int, int]],
) -> int:
    if not intervals:
        return 0

    intervals = sorted(intervals)

    total = 0
    start, end = intervals[0]

    for next_start, next_end in intervals[1:]:
        if next_start <= end + 1:
            end = max(end, next_end)
        else:
            total += end - start + 1
            start, end = next_start, next_end

    return total + end - start + 1


def build_summary(args) -> dict:
    genotypes = pl.read_csv(
        args.genotypes,
        separator="\t",
    )["genotype"].to_list()

    fasta_lengths = {
        path.stem: fasta_length(path)
        for path in args.clean_fastas
    }

    blocks = pl.read_csv(
        args.block_coords,
        separator="\t",
    ).with_columns(
        (
            pl.col("block_end_in_region")
            - pl.col("block_start_in_region")
            + 1
        ).alias("block_len_bp")
    )

    snps = pl.read_csv(
        args.snp_summary,
        separator="\t",
    )

    repeat_stats = pl.read_csv(
        args.masked_block_n_stats,
        separator="\t",
    )

    genotype_stats = {}

    for genotype in genotypes:
        genotype_blocks = blocks.filter(
            pl.col("genotype") == genotype
        )

        intervals = [
            (
                int(row["block_start_in_region"]),
                int(row["block_end_in_region"]),
            )
            for row in genotype_blocks.iter_rows(named=True)
        ]

        region_length = fasta_lengths[genotype]

        cumulated = (
            int(genotype_blocks["block_len_bp"].sum())
            if genotype_blocks.height
            else 0
        )

        covered = interval_union_length(intervals)

        genotype_repeats = repeat_stats.filter(
            pl.col("sample") == genotype
        )

        repeat_masked_bp = (
            int(
                genotype_repeats[
                    "repeat_masked_n_count"
                ].sum()
            )
            if genotype_repeats.height
            else 0
        )

        block_bp_for_repeat_stats = (
            int(
                genotype_repeats[
                    "unmasked_length_bp"
                ].sum()
            )
            if genotype_repeats.height
            else 0
        )

        repeat_masked_pct = (
            100
            * repeat_masked_bp
            / block_bp_for_repeat_stats
            if block_bp_for_repeat_stats
            else 0.0
        )

        genotype_stats[genotype] = {
            "region_length_bp": region_length,
            "cumulated_block_bp": cumulated,
            "covered_pct_of_region": round(
                100 * covered / region_length,
                2,
            ),
            "repeat_masked_pct_of_collinear_blocks": round(
                repeat_masked_pct,
                2,
            ),
        }

    diagnostic = snps.filter(
        pl.col("diagnostic_status") == "PASS"
    ).height

    summary = {
        "mode": args.mode,
        "input": {
            "n_genotypes": len(genotypes),
        },
        "repeat_masking": {
            "tool": "RepeatMasker",
            "simple_repeats_and_low_complexity": (
                args.repeat_masking
            ),
            "custom_library": (
                str(args.repeat_masking_library)
                if args.repeat_masking_library is not None
                else None
            ),
        },
        "global": {
            # Kept for compatibility with the current viewer.
            "n_blocks_kept": blocks["block_id"].n_unique(),
            "min_block_len_bp": int(
                blocks["block_len_bp"].min()
            ),
            "max_block_len_bp": int(
                blocks["block_len_bp"].max()
            ),
            "mean_block_len_bp": round(
                float(blocks["block_len_bp"].mean()),
                2,
            ),
            "n_snps_kept": snps.height,
        },
        "snp_discovery": {
            "detected_snps": snps.height,
            "diagnostic_snps": diagnostic,
            "non_diagnostic_snps": (
                snps.height - diagnostic
            ),
        },
        "failure_reasons": {
            "snps": count_reasons(
                snps,
                "final_status",
                "final_failure_reason",
            ),
        },
        "samples": genotype_stats,
    }

    if args.mode == "kasp":
        if args.assay_summary is None:
            raise ValueError(
                "--assay-summary is required in kasp mode"
            )

        assays = pl.read_csv(
            args.assay_summary,
            separator="\t",
        )

        passing_assays = assays.filter(
            pl.col("validation_status") == "PASS"
        ).height

        summary["kasp_assay_design"] = {
            "snps_with_assay_proposed": (
                assays["snp_id"].n_unique()
                if assays.height
                else 0
            ),
            "assays_proposed": assays.height,
        }

        summary["in_silico_validation"] = {
            "assays_passing_validation": (
                passing_assays
            ),
            "assays_failing_validation": (
                assays.height - passing_assays
            ),
        }

        summary["final_candidates"] = {
            "candidate_snps": snps.filter(
                pl.col("final_status") == "PASS"
            ).height,
            "candidate_assays": passing_assays,
        }

        summary["failure_reasons"]["assays"] = (
            count_reasons(
                assays,
                "validation_status",
                "validation_failure_reason",
            )
        )

    return summary


def add_section(
    lines: list[str],
    title: str,
):
    lines.extend([
        "",
        title,
        "-" * len(title),
    ])


def build_text(summary: dict) -> str:
    lines = [
        "KASPberry run summary",
        "=====================",
    ]

    add_section(lines, "Input")

    lines.append(
        "Genotypes (one input region each): "
        f"{summary['input']['n_genotypes']}"
    )

    masking = summary["repeat_masking"]

    add_section(lines, "Repeat masking")

    lines.append(
        "Simple repeats and low-complexity regions: "
        + (
            "masked"
            if masking[
                "simple_repeats_and_low_complexity"
            ]
            else "not masked"
        )
    )

    if masking["custom_library"]:
        lines.append(
            "Custom repeat library: "
            f"{masking['custom_library']}"
        )
    else:
        lines.append(
            "Custom repeat library: none"
        )

    stats = summary["global"]

    add_section(lines, "Collinear blocks")

    lines.extend([
        (
            "Collinear blocks: "
            f"{stats['n_blocks_kept']}"
        ),
        (
            "Smallest block length (bp): "
            f"{stats['min_block_len_bp']}"
        ),
        (
            "Largest block length (bp): "
            f"{stats['max_block_len_bp']}"
        ),
        (
            "Mean block length (bp): "
            f"{stats['mean_block_len_bp']:.2f}"
        ),
    ])

    stats = summary["snp_discovery"]

    add_section(lines, "SNP discovery")

    lines.extend([
        (
            "Detected SNPs: "
            f"{stats['detected_snps']}"
        ),
        (
            "Diagnostic SNPs: "
            f"{stats['diagnostic_snps']}"
        ),
        (
            "Non-diagnostic SNPs: "
            f"{stats['non_diagnostic_snps']}"
        ),
    ])

    if summary["mode"] == "kasp":
        stats = summary["kasp_assay_design"]

        add_section(
            lines,
            "KASP assay design",
        )

        lines.extend([
            (
                "SNPs with at least one assay proposed: "
                f"{stats['snps_with_assay_proposed']}"
            ),
            (
                "Assays proposed: "
                f"{stats['assays_proposed']}"
            ),
        ])

        stats = summary[
            "in_silico_validation"
        ]

        add_section(
            lines,
            "In silico validation",
        )

        lines.extend([
            (
                "Assays passing validation: "
                f"{stats['assays_passing_validation']}"
            ),
            (
                "Assays failing validation: "
                f"{stats['assays_failing_validation']}"
            ),
        ])

        stats = summary["final_candidates"]

        add_section(
            lines,
            "Final candidates",
        )

        lines.extend([
            (
                "Candidate SNPs: "
                f"{stats['candidate_snps']}"
            ),
            (
                "Candidate assays: "
                f"{stats['candidate_assays']}"
            ),
        ])

    reasons = summary["failure_reasons"]

    if any(reasons.values()):
        add_section(
            lines,
            "Failure reasons",
        )

        if reasons["snps"]:
            lines.append("SNPs:")

            lines.extend(
                f"  {reason}: {count}"
                for reason, count
                in reasons["snps"].items()
            )

        if reasons.get("assays"):
            lines.append("Assays:")

            lines.extend(
                f"  {reason}: {count}"
                for reason, count
                in reasons["assays"].items()
            )

    add_section(
        lines,
        "Per-genotype region statistics",
    )

    for genotype, stats in summary["samples"].items():
        lines.extend([
            genotype,
            (
                "  Region length (bp): "
                f"{stats['region_length_bp']}"
            ),
            (
                "  Cumulated collinear block "
                "length (bp): "
                f"{stats['cumulated_block_bp']}"
            ),
            (
                "  Region covered by collinear "
                "blocks (%): "
                f"{stats['covered_pct_of_region']:.2f}"
            ),
            (
                "  Repeat-masked bases in "
                "collinear blocks (%): "
                f"{stats['repeat_masked_pct_of_collinear_blocks']:.2f}"
            ),
            "",
        ])

    return "\n".join(lines).rstrip() + "\n"


def main():
    args = parse_args()

    summary = build_summary(args)

    args.json_output.parent.mkdir(
        parents=True,
        exist_ok=True,
    )
    args.txt_output.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    args.json_output.write_text(
        json.dumps(
            summary,
            indent=2,
        )
        + "\n"
    )

    args.txt_output.write_text(
        build_text(summary)
    )


if __name__ == "__main__":
    main()
