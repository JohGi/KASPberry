#!/usr/bin/env python3
"""Summarize MFEprimer in silico validation results for KASP assays."""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path


NOT_RUN = "NOT_RUN"

DIMER_HEADER = re.compile(r"^Dimer \d+:\s+(\S+)\s+x\s+(\S+)$")
DIMER_VALUES = re.compile(r"Score:\s+(\d+),\s+Delta G = (-?\d+(?:\.\d+)?)")
HAIRPIN_HEADER = re.compile(r"^Hairpin \d+:\s+(\S+)$")
HAIRPIN_VALUES = re.compile(
    r"Score:\s+(\d+),\s+Tm = (-?\d+(?:\.\d+)?) °C, "
    r"Delta G = (-?\d+(?:\.\d+)?)"
)
PRIMER_ASSAY = re.compile(
    r"^(snp::.+::assay::\d+)_(?:common|[ACGT]_specific)$"
)


@dataclass(frozen=True)
class ValidationContext:
    genotypes: list[str]
    design_by_snp: dict[str, dict[str, str]]
    design_by_snp_genotype: dict[tuple[str, str], dict[str, str]]
    positions_by_genotype: dict[tuple[str, str], tuple[str, int]]
    source_aliases: dict[str, str]
    canonical: dict[str, list[dict[str, str]]]
    noncanonical: dict[str, list[dict[str, str]]]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Summarize MFEprimer validation of PolyMarker KASP assays."
    )
    parser.add_argument(
        "--design-status",
        required=True,
        type=Path,
        help="PolyMarker SNP design status table.",
    )
    parser.add_argument(
        "--design-status-by-genotype",
        required=True,
        type=Path,
        help="PolyMarker SNP-by-genotype design status table.",
    )
    parser.add_argument(
        "--assays",
        required=True,
        type=Path,
        help="PolyMarker assay table.",
    )
    parser.add_argument(
        "--snp-positions",
        required=True,
        type=Path,
        help="Long SNP position table.",
    )
    parser.add_argument(
        "--aliases",
        nargs="+",
        required=True,
        type=Path,
        help="PolyMarker chromosome alias tables for KASP genotypes.",
    )
    parser.add_argument(
        "--in-silico-dir",
        required=True,
        type=Path,
        help="Directory containing MFEprimer specificity, dimer, and hairpin outputs.",
    )
    parser.add_argument(
        "--assay-status",
        required=True,
        type=Path,
        help="Output assay-level validation status table.",
    )
    parser.add_argument(
        "--assay-status-by-genotype",
        required=True,
        type=Path,
        help="Output assay-by-genotype specificity status table.",
    )
    parser.add_argument(
        "--validation-status",
        required=True,
        type=Path,
        help="Output SNP-level in silico validation status table.",
    )
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_source_aliases(paths: list[Path]) -> dict[str, str]:
    """Return the exact PolyMarker source-sequence alias for each genotype."""
    aliases: dict[str, str] = {}

    for path in paths:
        rows = read_tsv(path)
        source_rows = [
            row
            for row in rows
            if row["is_source_seq"].lower() == "true"
        ]

        if len(source_rows) != 1:
            raise ValueError(
                f"{path}: expected exactly one source sequence, "
                f"found {len(source_rows)}"
            )

        source = source_rows[0]
        genotype = source["genotype"]

        if genotype in aliases:
            raise ValueError(
                f"Duplicate chromosome alias table for genotype {genotype}"
            )

        aliases[genotype] = source["alias"]

    return aliases


def read_spec(path: Path) -> list[dict[str, str]]:
    """Read an MFEprimer .spec.tsv, including a header-only no-hit file."""
    with path.open(encoding="utf-8") as handle:
        lines = [
            line
            for line in handle
            if not line.startswith("#1-based coordinate")
        ]

    if not lines:
        return []

    if lines[0].startswith("#"):
        lines[0] = lines[0][1:]

    return list(csv.DictReader(lines, delimiter="\t"))


def write_tsv(
    path: Path,
    fields: list[str],
    rows: list[dict[str, object]],
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


def primer_pair_id(primer_name: str) -> str:
    """Return pair ID from an MFEprimer fp/rp primer name."""
    if primer_name.endswith("_fp"):
        return primer_name.removesuffix("_fp")

    if primer_name.endswith("_rp"):
        return primer_name.removesuffix("_rp")

    raise ValueError(
        f"Unexpected MFEprimer primer name: {primer_name}"
    )


def pair_id(hit: dict[str, str]) -> str:
    fp = primer_pair_id(hit["fpName"])
    rp = primer_pair_id(hit["rpName"])

    if fp != rp:
        raise ValueError(
            "Inconsistent primer-pair names in MFEprimer output: "
            f"{hit['fpName']} / {hit['rpName']}"
        )

    return fp


def assay_from_primer(primer_id: str) -> str | None:
    match = PRIMER_ASSAY.match(primer_id)
    return match.group(1) if match else None


def parse_dimers(path: Path) -> list[tuple[str, str, int, float]]:
    records: list[tuple[str, str, int, float]] = []
    pending: tuple[str, str] | None = None

    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()

            match = DIMER_HEADER.match(line)
            if match:
                pending = (match.group(1), match.group(2))
                continue

            match = DIMER_VALUES.search(line)
            if match and pending:
                records.append(
                    (
                        pending[0],
                        pending[1],
                        int(match.group(1)),
                        float(match.group(2)),
                    )
                )
                pending = None

    return records


def parse_hairpins(path: Path) -> list[tuple[str, int, float, float]]:
    records: list[tuple[str, int, float, float]] = []
    pending: str | None = None

    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()

            match = HAIRPIN_HEADER.match(line)
            if match:
                pending = match.group(1)
                continue

            match = HAIRPIN_VALUES.search(line)
            if match and pending:
                records.append(
                    (
                        pending,
                        int(match.group(1)),
                        float(match.group(2)),
                        float(match.group(3)),
                    )
                )
                pending = None

    return records


def target_hit(
    hit: dict[str, str],
    target_alias: str,
    target_position: int,
) -> bool:
    """Check that the allele-specific primer ends exactly on the target SNP."""
    if hit["chrom"] != target_alias:
        return False

    fp_name = hit["fpName"]
    rp_name = hit["rpName"]
    fp_is_input_primer = fp_name.endswith("_fp")
    rp_is_input_primer = rp_name.endswith("_fp")
    fp_is_pair_partner = fp_name.endswith("_rp")
    rp_is_pair_partner = rp_name.endswith("_rp")

    if (
        fp_is_input_primer + rp_is_input_primer != 1
        or fp_is_pair_partner + rp_is_pair_partner != 1
    ):
        raise ValueError(
            "Expected exactly one _fp and one _rp member in MFEprimer "
            f"hit: {fp_name} / {rp_name}"
        )

    # The first primer supplied in each canonical pair is allele-specific.
    # MFEprimer keeps that input role as the "_fp" suffix even when genomic
    # orientation places it in the rpName/rpStart output columns.
    if fp_is_input_primer:
        return int(hit["fpEnd"]) == target_position

    return int(hit["rpStart"]) == target_position


def build_validation_context(
    args: argparse.Namespace,
    design_status: list[dict[str, str]],
) -> ValidationContext:
    """Load and index the per-genotype inputs needed for validation."""
    by_genotype = read_tsv(args.design_status_by_genotype)
    positions = read_tsv(args.snp_positions)

    genotypes = sorted({
        row["genotype"]
        for row in by_genotype
    })

    source_aliases = read_source_aliases(args.aliases)
    if set(source_aliases) != set(genotypes):
        raise ValueError(
            "Genotypes differ between PolyMarker alias tables and "
            "in silico validation inputs"
        )

    return ValidationContext(
        genotypes=genotypes,
        design_by_snp={
            row["snp_id"]: row
            for row in design_status
        },
        design_by_snp_genotype={
            (row["snp_id"], row["genotype"]): row
            for row in by_genotype
        },
        positions_by_genotype={
            (row["snp_id"], row["genotype"]): (
                row["nt"],
                int(row["pos_in_source_seq"]),
            )
            for row in positions
        },
        source_aliases=source_aliases,
        canonical={
            genotype: read_spec(
                args.in_silico_dir
                / "specificity"
                / genotype
                / "canonical.spec.tsv"
            )
            for genotype in genotypes
        },
        noncanonical={
            genotype: read_spec(
                args.in_silico_dir
                / "specificity"
                / genotype
                / "noncanonical.spec.tsv"
            )
            for genotype in genotypes
        },
    )


def find_bad_dimer_assays(path: Path) -> set[str]:
    """Return assays with at least one reported intra-assay dimer."""
    bad_assays: set[str] = set()

    for primer_a, primer_b, _score, _dg in parse_dimers(path):
        assay_a = assay_from_primer(primer_a)
        assay_b = assay_from_primer(primer_b)

        # Cross-assay dimers are intentionally ignored.
        if assay_a and assay_a == assay_b:
            bad_assays.add(assay_a)

    return bad_assays


def find_bad_hairpin_assays(path: Path) -> set[str]:
    """Return assays with at least one reported primer hairpin."""
    bad_assays: set[str] = set()

    for primer, _score, _tm, _dg in parse_hairpins(path):
        assay_id = assay_from_primer(primer)
        if assay_id:
            bad_assays.add(assay_id)

    return bad_assays


def evaluate_assay_in_genotype(
    assay: dict[str, str],
    genotype: str,
    context: ValidationContext,
) -> tuple[dict[str, object], list[str]]:
    """Evaluate specificity of one assay in one genotype."""
    assay_id = assay["assay_id"]
    snp_id = assay["snp_id"]
    design_key = (snp_id, genotype)

    if design_key not in context.design_by_snp_genotype:
        raise ValueError(
            f"{snp_id}/{genotype} is missing from the "
            "PolyMarker by-genotype status table"
        )

    allele = context.design_by_snp_genotype[design_key]["expected_allele"]

    position_key = (snp_id, genotype)
    if position_key not in context.positions_by_genotype:
        raise ValueError(
            f"{snp_id}/{genotype} is missing from snp_positions_long.tsv"
        )

    position_allele, position = context.positions_by_genotype[position_key]
    if allele != position_allele:
        raise ValueError(
            f"{snp_id}/{genotype}: expected allele {allele} != "
            f"snp_positions_long.tsv allele {position_allele}"
        )

    allele1 = assay["first_allele"]
    allele2 = assay["second_allele"]
    if allele == allele1:
        genotype_allele_id = "1"
    elif allele == allele2:
        genotype_allele_id = "2"
    else:
        raise ValueError(
            f"{snp_id}/{genotype}: allele {allele} is absent from {assay_id}"
        )

    pair_definitions = {
        "allele1/common_amplicons": (
            f"{assay_id}_{allele1}_common",
            context.canonical[genotype],
        ),
        "allele2/common_amplicons": (
            f"{assay_id}_{allele2}_common",
            context.canonical[genotype],
        ),
        "allele1/allele2_amplicons": (
            f"{assay_id}_{allele1}_{allele2}_allele_pair",
            context.noncanonical[genotype],
        ),
        "allele1/allele1_amplicons": (
            f"{assay_id}_{allele1}_{allele1}_allele_self",
            context.noncanonical[genotype],
        ),
        "allele2/allele2_amplicons": (
            f"{assay_id}_{allele2}_{allele2}_allele_self",
            context.noncanonical[genotype],
        ),
        "common/common_amplicons": (
            f"{assay_id}_common_self",
            context.noncanonical[genotype],
        ),
    }
    pair_counts = {
        column: sum(
            pair_id(hit) == expected_pair_id
            for hit in hits
        )
        for column, (expected_pair_id, hits) in pair_definitions.items()
    }
    expected_pair = f"{assay_id}_{allele}_common"
    expected_hits = [
        hit
        for hit in context.canonical[genotype]
        if pair_id(hit) == expected_pair
    ]
    target_hits = [
        hit
        for hit in expected_hits
        if target_hit(
            hit,
            context.source_aliases[genotype],
            position,
        )
    ]

    failure_reasons: list[str] = []
    if not target_hits:
        failure_reasons.append("missing_target_amplicon")
    total_amplicons = sum(pair_counts.values())
    unexpected_count = total_amplicons - min(len(target_hits), 1)
    if unexpected_count > 0:
        failure_reasons.append("unexpected_amplicon")

    return (
        {
            "assay_id": assay_id,
            "genotype": genotype,
            "allele1": allele1,
            "allele2": allele2,
            "genotype_allele_id": genotype_allele_id,
            "target_amplicons": len(target_hits),
            **pair_counts,
            "status": "PASS" if not failure_reasons else "FAIL",
            "failure_reason": ";".join(failure_reasons),
        },
        failure_reasons,
    )


def evaluate_assays(
    assays: list[dict[str, str]],
    context: ValidationContext,
    bad_dimer_assays: set[str],
    bad_hairpin_assays: set[str],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """Evaluate all assays and return assay- and genotype-level rows."""
    assay_rows: list[dict[str, object]] = []
    genotype_rows: list[dict[str, object]] = []

    for assay in assays:
        assay_id = assay["assay_id"]
        snp_id = assay["snp_id"]

        if snp_id not in context.design_by_snp:
            raise ValueError(
                f"{assay_id}: {snp_id} is missing from the "
                "PolyMarker design status table"
            )

        failure_reasons: set[str] = set()

        for genotype in context.genotypes:
            genotype_row, genotype_failures = evaluate_assay_in_genotype(
                assay,
                genotype,
                context,
            )
            genotype_rows.append(genotype_row)
            failure_reasons.update(genotype_failures)

        if assay_id in bad_dimer_assays:
            failure_reasons.add("dimer")
        if assay_id in bad_hairpin_assays:
            failure_reasons.add("hairpin")

        assay_rows.append({
            "assay_id": assay_id,
            "snp_id": snp_id,
            "status": "PASS" if not failure_reasons else "FAIL",
            "failure_reason": ";".join(sorted(failure_reasons)),
        })

    return assay_rows, genotype_rows


def build_snp_validation_rows(
    design_status: list[dict[str, str]],
    assays: list[dict[str, str]],
    assay_rows: list[dict[str, object]],
) -> list[dict[str, object]]:
    """Roll assay-level validation results up to one status per SNP."""
    assay_status_by_id = {
        row["assay_id"]: row
        for row in assay_rows
    }

    assays_by_snp: dict[str, list[str]] = defaultdict(list)
    for assay in assays:
        assays_by_snp[assay["snp_id"]].append(assay["assay_id"])

    snp_rows: list[dict[str, object]] = []

    for design_row in design_status:
        snp_id = design_row["snp_id"]

        if design_row["status"] != "PASS":
            status = NOT_RUN
            failure_reason = ""
        else:
            assay_ids = assays_by_snp.get(snp_id, [])
            if not assay_ids:
                raise ValueError(
                    f"{snp_id} passed PolyMarker design but has no assay "
                    "in the PolyMarker assay table"
                )

            if any(
                assay_status_by_id[assay_id]["status"] == "PASS"
                for assay_id in assay_ids
            ):
                status = "PASS"
                failure_reason = ""
            else:
                status = "FAIL"
                failure_reason = "no_assay_passed_in_silico_validation"

        snp_rows.append({
            "snp_id": snp_id,
            "status": status,
            "failure_reason": failure_reason,
        })

    return snp_rows


def main() -> None:
    args = parse_args()

    design_status = read_tsv(args.design_status)
    assays = read_tsv(args.assays)
    context = build_validation_context(args, design_status)

    bad_dimer_assays = find_bad_dimer_assays(
        args.in_silico_dir / "dimers.tsv"
    )
    bad_hairpin_assays = find_bad_hairpin_assays(
        args.in_silico_dir / "hairpins.tsv"
    )

    assay_rows, genotype_rows = evaluate_assays(
        assays,
        context,
        bad_dimer_assays,
        bad_hairpin_assays,
    )
    snp_rows = build_snp_validation_rows(
        design_status,
        assays,
        assay_rows,
    )

    write_tsv(
        args.assay_status,
        ["assay_id", "snp_id", "status", "failure_reason"],
        assay_rows,
    )
    write_tsv(
        args.assay_status_by_genotype,
        [
            "assay_id",
            "genotype",
            "allele1",
            "allele2",
            "genotype_allele_id",
            "target_amplicons",
            "allele1/common_amplicons",
            "allele2/common_amplicons",
            "allele1/allele2_amplicons",
            "allele1/allele1_amplicons",
            "allele2/allele2_amplicons",
            "common/common_amplicons",
            "status",
            "failure_reason",
        ],
        genotype_rows,
    )
    write_tsv(
        args.validation_status,
        ["snp_id", "status", "failure_reason"],
        snp_rows,
    )


if __name__ == "__main__":
    main()
