#!/usr/bin/env python3
# Author: Johanna Girodolle

"""Read aggregated SNP and KASP assay results for the region viewer."""

from __future__ import annotations

import csv
from pathlib import Path

from .models import AssayResult, SnpResult


COMMON_SNP_RESULT_COLUMNS = {
    "snp_id",
    "block_id",
    "aln_pos",
    "diagnostic_status",
    "diagnostic_failure_reason",
    "final_status",
    "final_failure_reason",
}

KASP_SNP_RESULT_COLUMNS = {
    "design_status",
    "design_failure_reason",
    "validation_status",
    "validation_failure_reason",
}

ASSAY_RESULT_COLUMNS = {
    "assay_id",
    "snp_id",
    "first_allele",
    "second_allele",
    "first_primer",
    "second_primer",
    "common_primer",
    "first_primer_with_tail",
    "second_primer_with_tail",
    "source_genotypes",
    "validation_status",
    "validation_failure_reason",
}


def read_tsv_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    """Read a TSV file and return its header and rows."""
    if not path.is_file():
        raise FileNotFoundError(f"TSV file not found: {path}")

    with path.open(
        newline="",
        encoding="utf-8",
    ) as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        if reader.fieldnames is None:
            raise ValueError(f"TSV file has no header: {path}")

        return list(reader.fieldnames), list(reader)


def require_columns(
    path: Path,
    fieldnames: list[str],
    required_columns: set[str],
) -> None:
    """Ensure that all required columns are present."""
    missing_columns = required_columns - set(fieldnames)

    if missing_columns:
        raise ValueError(
            f"Missing required columns in {path}: "
            + ", ".join(sorted(missing_columns))
        )


def read_snp_results(
    path: Path,
    mode: str,
) -> dict[str, SnpResult]:
    """Read aggregated workflow status for every detected SNP."""
    if mode not in {"snps", "kasp"}:
        raise ValueError(f"Unsupported viewer mode: {mode}")

    fieldnames, rows = read_tsv_rows(path)

    required_columns = set(COMMON_SNP_RESULT_COLUMNS)

    if mode == "kasp":
        required_columns.update(KASP_SNP_RESULT_COLUMNS)

    require_columns(
        path=path,
        fieldnames=fieldnames,
        required_columns=required_columns,
    )

    results: dict[str, SnpResult] = {}

    for row in rows:
        snp_id = row["snp_id"]

        if snp_id in results:
            raise ValueError(
                f"Duplicate snp_id in {path}: {snp_id}"
            )

        result = SnpResult(
            snp_id=snp_id,
            block_id=row["block_id"],
            aln_pos=int(row["aln_pos"]),
            diagnostic_status=row["diagnostic_status"],
            diagnostic_failure_reason=row[
                "diagnostic_failure_reason"
            ],
            final_status=row["final_status"],
            final_failure_reason=row["final_failure_reason"],
            design_status=(
                row["design_status"]
                if mode == "kasp"
                else None
            ),
            design_failure_reason=(
                row["design_failure_reason"]
                if mode == "kasp"
                else None
            ),
            validation_status=(
                row["validation_status"]
                if mode == "kasp"
                else None
            ),
            validation_failure_reason=(
                row["validation_failure_reason"]
                if mode == "kasp"
                else None
            ),
        )

        results[snp_id] = result

    return results


def read_assay_results(
    path: Path,
) -> dict[str, AssayResult]:
    """Read all PolyMarker assays and their in silico validation status."""
    fieldnames, rows = read_tsv_rows(path)

    require_columns(
        path=path,
        fieldnames=fieldnames,
        required_columns=ASSAY_RESULT_COLUMNS,
    )

    results: dict[str, AssayResult] = {}

    for row in rows:
        assay_id = row["assay_id"]

        if assay_id in results:
            raise ValueError(
                f"Duplicate assay_id in {path}: {assay_id}"
            )

        results[assay_id] = AssayResult(
            assay_id=assay_id,
            snp_id=row["snp_id"],
            first_allele=row["first_allele"],
            second_allele=row["second_allele"],
            first_primer=row["first_primer"],
            second_primer=row["second_primer"],
            common_primer=row["common_primer"],
            first_primer_with_tail=row[
                "first_primer_with_tail"
            ],
            second_primer_with_tail=row[
                "second_primer_with_tail"
            ],
            source_genotypes=row["source_genotypes"],
            validation_status=row["validation_status"],
            validation_failure_reason=row[
                "validation_failure_reason"
            ],
        )

    return results


def group_assays_by_snp(
    assays: dict[str, AssayResult],
) -> dict[str, list[AssayResult]]:
    """Group assay results by their targeted SNP."""
    assays_by_snp: dict[str, list[AssayResult]] = {}

    for assay in assays.values():
        assays_by_snp.setdefault(
            assay.snp_id,
            [],
        ).append(assay)

    for snp_assays in assays_by_snp.values():
        snp_assays.sort(
            key=lambda assay: assay.assay_id
        )

    return assays_by_snp
