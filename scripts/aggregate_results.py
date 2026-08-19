#!/usr/bin/env python3
"""Aggregate KASPberry SNP and assay results into user-facing output files."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from ids import make_snp_id


NOT_RUN = "NOT_RUN"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Aggregate KASPberry results for the SNP or KASP workflow mode."
    )
    parser.add_argument(
        "--mode",
        required=True,
        choices=["snps", "kasp"],
        help="Workflow mode to summarize.",
    )
    parser.add_argument(
        "--snps-vcf",
        required=True,
        type=Path,
        help="VCF containing all detected SNPs.",
    )
    parser.add_argument(
        "--snp-positions",
        required=True,
        type=Path,
        help="Wide SNP position table containing all detected SNPs.",
    )
    parser.add_argument(
        "--diagnostic-status",
        required=True,
        type=Path,
        help="Diagnostic SNP status table.",
    )
    parser.add_argument(
        "--design-status",
        type=Path,
        help="PolyMarker design status table. Required in kasp mode.",
    )
    parser.add_argument(
        "--validation-status",
        type=Path,
        help="In silico SNP validation status table. Required in kasp mode.",
    )
    parser.add_argument(
        "--assays",
        type=Path,
        help="PolyMarker assay table. Required in kasp mode.",
    )
    parser.add_argument(
        "--assay-status",
        type=Path,
        help="In silico assay status table. Required in kasp mode.",
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Output SNP summary TSV.",
    )
    parser.add_argument(
        "--retained-vcf",
        required=True,
        type=Path,
        help=(
            "Output VCF containing retained SNPs: diagnostic SNPs in snps mode, "
            "candidate SNPs in kasp mode."
        ),
    )
    parser.add_argument(
        "--assay-summary",
        type=Path,
        help="Output assay summary TSV. Required in kasp mode.",
    )
    parser.add_argument(
        "--candidate-assays",
        type=Path,
        help="Output final candidate assay TSV. Required in kasp mode.",
    )
    parser.add_argument(
        "--primers-to-order",
        type=Path,
        help="Output primer-ordering TSV. Required in kasp mode.",
    )
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def index_by(
    rows: list[dict[str, str]],
    key: str,
) -> dict[str, dict[str, str]]:
    indexed: dict[str, dict[str, str]] = {}

    for row in rows:
        value = row[key]
        if value in indexed:
            raise ValueError(f"Duplicate {key}: {value}")
        indexed[value] = row

    return indexed


def write_tsv(
    path: Path,
    fields: list[str],
    rows: list[dict[str, str]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def write_filtered_vcf(
    input_vcf: Path,
    output_vcf: Path,
    retained_snp_ids: set[str],
) -> None:
    """Copy VCF headers and retain only records whose canonical SNP ID is selected."""
    output_vcf.parent.mkdir(parents=True, exist_ok=True)

    found: set[str] = set()

    with (
        input_vcf.open(encoding="utf-8") as fin,
        output_vcf.open("w", encoding="utf-8") as fout,
    ):
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                raise ValueError(f"Malformed VCF line in {input_vcf}: {line.rstrip()}")

            block_id = fields[0].removesuffix(".aln")
            snp_id = make_snp_id(block_id, fields[1])

            if snp_id in retained_snp_ids:
                fout.write(line)
                found.add(snp_id)

    missing = retained_snp_ids - found
    if missing:
        preview = ", ".join(sorted(missing)[:10])
        suffix = " ..." if len(missing) > 10 else ""
        raise ValueError(
            f"{len(missing)} retained SNP(s) were not found in {input_vcf}: "
            f"{preview}{suffix}"
        )


def require_row(
    table: dict[str, dict[str, str]],
    key: str,
    table_name: str,
) -> dict[str, str]:
    try:
        return table[key]
    except KeyError as exc:
        raise ValueError(
            f"{key} is present in SNP positions but missing from {table_name}"
        ) from exc


def aggregate_snps(
    positions: list[dict[str, str]],
    diagnostic_by_snp: dict[str, dict[str, str]],
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []

    for position_row in positions:
        snp_id = position_row["snp_id"]
        diagnostic = require_row(
            diagnostic_by_snp,
            snp_id,
            "diagnostic status",
        )

        diagnostic_status = diagnostic["status"]
        diagnostic_failure_reason = diagnostic["failure_reason"]

        rows.append({
            **position_row,
            "diagnostic_status": diagnostic_status,
            "diagnostic_failure_reason": diagnostic_failure_reason,
            "final_status": diagnostic_status,
            "final_failure_reason": diagnostic_failure_reason,
        })

    return rows


def aggregate_kasp(
    positions: list[dict[str, str]],
    diagnostic_by_snp: dict[str, dict[str, str]],
    design_by_snp: dict[str, dict[str, str]],
    validation_by_snp: dict[str, dict[str, str]],
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []

    for position_row in positions:
        snp_id = position_row["snp_id"]
        diagnostic = require_row(
            diagnostic_by_snp,
            snp_id,
            "diagnostic status",
        )

        diagnostic_status = diagnostic["status"]
        diagnostic_failure_reason = diagnostic["failure_reason"]

        if diagnostic_status != "PASS":
            design_status = NOT_RUN
            design_failure_reason = ""
            validation_status = NOT_RUN
            validation_failure_reason = ""
            final_status = "FAIL"
            final_failure_reason = diagnostic_failure_reason

        else:
            design = design_by_snp.get(snp_id)

            if design is None:
                raise ValueError(
                    f"{snp_id} passed diagnostic filtering but is missing "
                    "from the PolyMarker design status table"
                )

            design_status = design["status"]
            design_failure_reason = design["failure_reason"]

            if design_status != "PASS":
                validation_status = NOT_RUN
                validation_failure_reason = ""
                final_status = "FAIL"
                final_failure_reason = design_failure_reason

            else:
                validation = validation_by_snp.get(snp_id)

                if validation is None:
                    raise ValueError(
                        f"{snp_id} passed PolyMarker design but is missing "
                        "from the in silico validation status table"
                    )

                validation_status = validation["status"]
                validation_failure_reason = validation["failure_reason"]
                final_status = validation_status
                final_failure_reason = validation_failure_reason

        rows.append({
            **position_row,
            "diagnostic_status": diagnostic_status,
            "diagnostic_failure_reason": diagnostic_failure_reason,
            "design_status": design_status,
            "design_failure_reason": design_failure_reason,
            "validation_status": validation_status,
            "validation_failure_reason": validation_failure_reason,
            "final_status": final_status,
            "final_failure_reason": final_failure_reason,
        })

    return rows


def aggregate_assays(
    assays: list[dict[str, str]],
    assay_status_by_id: dict[str, dict[str, str]],
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []

    for assay in assays:
        assay_id = assay["assay_id"]

        try:
            status = assay_status_by_id[assay_id]
        except KeyError as exc:
            raise ValueError(
                f"{assay_id} is present in the PolyMarker assay table but "
                "missing from the in silico assay status table"
            ) from exc

        rows.append({
            **assay,
            "validation_status": status["status"],
            "validation_failure_reason": status["failure_reason"],
        })

    return rows


def build_primers_to_order(
    candidate_assays: list[dict[str, str]],
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []

    for assay in candidate_assays:
        assay_id = assay["assay_id"]
        snp_id = assay["snp_id"]

        rows.extend([
            {
                "assay_id": assay_id,
                "snp_id": snp_id,
                "primer_role": "first",
                "allele": assay["first_allele"],
                "sequence": assay["first_primer_with_tail"],
            },
            {
                "assay_id": assay_id,
                "snp_id": snp_id,
                "primer_role": "second",
                "allele": assay["second_allele"],
                "sequence": assay["second_primer_with_tail"],
            },
            {
                "assay_id": assay_id,
                "snp_id": snp_id,
                "primer_role": "common",
                "allele": "",
                "sequence": assay["common_primer"],
            },
        ])

    return rows


def require_kasp_args(args: argparse.Namespace) -> None:
    required = {
        "--design-status": args.design_status,
        "--validation-status": args.validation_status,
        "--assays": args.assays,
        "--assay-status": args.assay_status,
        "--assay-summary": args.assay_summary,
        "--candidate-assays": args.candidate_assays,
        "--primers-to-order": args.primers_to_order,
    }

    missing = [
        option
        for option, value in required.items()
        if value is None
    ]

    if missing:
        raise ValueError(
            "The following arguments are required in kasp mode: "
            + ", ".join(missing)
        )


def main() -> None:
    args = parse_args()

    positions = read_tsv(args.snp_positions)
    diagnostic_by_snp = index_by(
        read_tsv(args.diagnostic_status),
        "snp_id",
    )

    position_fields = list(positions[0]) if positions else [
        "snp_id",
        "block_id",
        "aln_pos",
    ]

    if args.mode == "snps":
        snp_rows = aggregate_snps(
            positions=positions,
            diagnostic_by_snp=diagnostic_by_snp,
        )

        fields = [
            *position_fields,
            "diagnostic_status",
            "diagnostic_failure_reason",
            "final_status",
            "final_failure_reason",
        ]

        write_tsv(args.output, fields, snp_rows)

        retained_snp_ids = {
            row["snp_id"]
            for row in snp_rows
            if row["final_status"] == "PASS"
        }
        write_filtered_vcf(
            args.snps_vcf,
            args.retained_vcf,
            retained_snp_ids,
        )
        return

    require_kasp_args(args)

    design_by_snp = index_by(
        read_tsv(args.design_status),
        "snp_id",
    )
    validation_by_snp = index_by(
        read_tsv(args.validation_status),
        "snp_id",
    )

    snp_rows = aggregate_kasp(
        positions=positions,
        diagnostic_by_snp=diagnostic_by_snp,
        design_by_snp=design_by_snp,
        validation_by_snp=validation_by_snp,
    )

    snp_fields = [
        *position_fields,
        "diagnostic_status",
        "diagnostic_failure_reason",
        "design_status",
        "design_failure_reason",
        "validation_status",
        "validation_failure_reason",
        "final_status",
        "final_failure_reason",
    ]

    write_tsv(args.output, snp_fields, snp_rows)

    retained_snp_ids = {
        row["snp_id"]
        for row in snp_rows
        if row["final_status"] == "PASS"
    }
    write_filtered_vcf(
        args.snps_vcf,
        args.retained_vcf,
        retained_snp_ids,
    )

    assays = read_tsv(args.assays)
    assay_status_by_id = index_by(
        read_tsv(args.assay_status),
        "assay_id",
    )

    assay_rows = aggregate_assays(
        assays=assays,
        assay_status_by_id=assay_status_by_id,
    )

    assay_fields = (
        [
            *list(assays[0]),
            "validation_status",
            "validation_failure_reason",
        ]
        if assays
        else [
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
        ]
    )

    write_tsv(
        args.assay_summary,
        assay_fields,
        assay_rows,
    )

    candidate_assays = [
        row
        for row in assay_rows
        if row["validation_status"] == "PASS"
    ]

    write_tsv(
        args.candidate_assays,
        assay_fields,
        candidate_assays,
    )

    primers_to_order = build_primers_to_order(candidate_assays)

    write_tsv(
        args.primers_to_order,
        [
            "assay_id",
            "snp_id",
            "primer_role",
            "allele",
            "sequence",
        ],
        primers_to_order,
    )


if __name__ == "__main__":
    main()
