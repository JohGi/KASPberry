#!/usr/bin/env python3
"""Summarize MFEprimer in silico validation results for KASP assays."""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path


NOT_RUN = "NOT_RUN"

DIMER_HEADER = re.compile(r"^Dimer \d+:\s+(\S+)\s+x\s+(\S+)$")
DIMER_VALUES = re.compile(r"Score:\s+(\d+),\s+Delta G = (-?\d+(?:\.\d+)?)")
HAIRPIN_HEADER = re.compile(r"^Hairpin \d+:\s+(\S+)$")
HAIRPIN_VALUES = re.compile(
    r"Score:\s+(\d+),\s+Tm = (-?\d+(?:\.\d+)?) °C, "
    r"Delta G = (-?\d+(?:\.\d+)?)"
)
PRIMER_ASSAY = re.compile(r"^(snp::.+::assay::\d+)_(?:comm_rev|[ACGT]_fw)$")


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


def pair_id(hit: dict[str, str]) -> str:
    fp = hit["fpName"].removesuffix("_fp")
    rp = hit["rpName"].removesuffix("_rp")

    if fp != rp:
        raise ValueError(
            "Inconsistent primer-pair names in MFEprimer output: "
            f"{fp} / {rp}"
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
    target_chromosome: str,
    target_position: int,
    allele_primer: str,
) -> bool:
    """Check that the allele-specific primer ends exactly on the target SNP."""
    if hit["chrom"].split("__", 1)[0] != target_chromosome:
        return False

    if hit["fpSeq"].upper() == allele_primer:
        return int(hit["fpEnd"]) == target_position

    if hit["rpSeq"].upper() == allele_primer:
        return int(hit["rpEnd"]) == target_position

    return False


def main() -> None:
    args = parse_args()

    design_status = read_tsv(args.design_status)
    by_genotype = read_tsv(args.design_status_by_genotype)
    assays = read_tsv(args.assays)
    positions = read_tsv(args.snp_positions)

    design_by_snp = {
        row["snp_id"]: row
        for row in design_status
    }

    design_by_snp_genotype = {
        (row["snp_id"], row["genotype"]): row
        for row in by_genotype
    }

    positions_by_genotype = {
        (row["snp_id"], row["genotype"]): (
            row["nt"],
            int(row["pos_in_source_seq"]),
        )
        for row in positions
    }

    genotypes = sorted({
        row["genotype"]
        for row in by_genotype
    })

    canonical = {
        genotype: read_spec(
            args.in_silico_dir
            / "specificity"
            / genotype
            / "canonical.spec.tsv"
        )
        for genotype in genotypes
    }

    noncanonical = {
        genotype: read_spec(
            args.in_silico_dir
            / "specificity"
            / genotype
            / "noncanonical.spec.tsv"
        )
        for genotype in genotypes
    }

    bad_dimers: dict[
        str,
        list[tuple[str, str, int, float]],
    ] = defaultdict(list)

    for primer_a, primer_b, score, dg in parse_dimers(
        args.in_silico_dir / "dimers.tsv"
    ):
        assay_a = assay_from_primer(primer_a)
        assay_b = assay_from_primer(primer_b)

        # Cross-assay dimers are intentionally ignored.
        if assay_a and assay_a == assay_b:
            bad_dimers[assay_a].append(
                (primer_a, primer_b, score, dg)
            )

    bad_hairpins: dict[
        str,
        list[tuple[str, int, float, float]],
    ] = defaultdict(list)

    for primer, score, tm, dg in parse_hairpins(
        args.in_silico_dir / "hairpins.tsv"
    ):
        assay_id = assay_from_primer(primer)

        if assay_id:
            bad_hairpins[assay_id].append(
                (primer, score, tm, dg)
            )

    assay_rows: list[dict[str, object]] = []
    genotype_rows: list[dict[str, object]] = []

    for assay in assays:
        assay_id = assay["assay_id"]
        snp_id = assay["snp_id"]

        if snp_id not in design_by_snp:
            raise ValueError(
                f"{assay_id}: {snp_id} is missing from the "
                "PolyMarker design status table"
            )

        specificity_results: list[bool] = []
        assay_failure_reasons: set[str] = set()

        for genotype in genotypes:
            design_key = (snp_id, genotype)

            if design_key not in design_by_snp_genotype:
                raise ValueError(
                    f"{snp_id}/{genotype} is missing from the "
                    "PolyMarker by-genotype status table"
                )

            design_row = design_by_snp_genotype[design_key]
            allele = design_row["expected_allele"]
            target_chromosome = design_row["target_chromosome"]

            position_key = (snp_id, genotype)
            if position_key not in positions_by_genotype:
                raise ValueError(
                    f"{snp_id}/{genotype} is missing from "
                    "snp_positions_long.tsv"
                )

            position_allele, position = positions_by_genotype[position_key]

            if allele != position_allele:
                raise ValueError(
                    f"{snp_id}/{genotype}: expected allele {allele} != "
                    f"snp_positions_long.tsv allele {position_allele}"
                )

            if allele == assay["first_allele"]:
                allele_primer = assay["first_primer"].upper()
                other_allele = assay["second_allele"]

            elif allele == assay["second_allele"]:
                allele_primer = assay["second_primer"].upper()
                other_allele = assay["first_allele"]

            else:
                raise ValueError(
                    f"{snp_id}/{genotype}: allele {allele} "
                    f"is absent from {assay_id}"
                )

            expected_pair = f"{assay_id}_{allele}_common"
            unexpected_pair = f"{assay_id}_{other_allele}_common"

            expected_hits = [
                hit
                for hit in canonical[genotype]
                if pair_id(hit) == expected_pair
            ]

            unexpected_hits = [
                hit
                for hit in canonical[genotype]
                if pair_id(hit) == unexpected_pair
            ]

            noncanonical_hits = [
                hit
                for hit in noncanonical[genotype]
                if pair_id(hit).startswith(f"{assay_id}_")
            ]

            target_hits = [
                hit
                for hit in expected_hits
                if target_hit(
                    hit,
                    target_chromosome,
                    position,
                    allele_primer,
                )
            ]

            failure_reasons: list[str] = []

            if len(target_hits) == 0:
                failure_reasons.append("missing_target_amplicon")

            if len(expected_hits) > 1:
                failure_reasons.append("multiple_expected_amplicons")

            if unexpected_hits:
                failure_reasons.append("unexpected_allele_amplicon")

            if noncanonical_hits:
                failure_reasons.append("noncanonical_amplicon")

            passed = not failure_reasons

            specificity_results.append(passed)
            assay_failure_reasons.update(failure_reasons)

            genotype_rows.append({
                "assay_id": assay_id,
                "genotype": genotype,
                "expected_allele": allele,
                "status": "PASS" if passed else "FAIL",
                "failure_reason": ";".join(failure_reasons),
                "expected_amplicons": len(expected_hits),
                "target_amplicons": len(target_hits),
                "unexpected_amplicons": len(unexpected_hits),
                "noncanonical_amplicons": len(noncanonical_hits),
            })

        if bad_dimers[assay_id]:
            assay_failure_reasons.add("dimer")

        if bad_hairpins[assay_id]:
            assay_failure_reasons.add("hairpin")

        validation_pass = not assay_failure_reasons

        assay_rows.append({
            "assay_id": assay_id,
            "snp_id": snp_id,
            "status": "PASS" if validation_pass else "FAIL",
            "failure_reason": ";".join(
                sorted(assay_failure_reasons)
            ),
        })

    assay_status_by_id = {
        row["assay_id"]: row
        for row in assay_rows
    }

    assays_by_snp: dict[str, list[str]] = defaultdict(list)
    for assay in assays:
        assays_by_snp[assay["snp_id"]].append(
            assay["assay_id"]
        )

    snp_rows: list[dict[str, object]] = []

    for row in design_status:
        snp_id = row["snp_id"]

        if row["status"] != "PASS":
            status = NOT_RUN
            failure_reason = ""

        else:
            assay_ids = assays_by_snp.get(snp_id, [])

            if not assay_ids:
                raise ValueError(
                    f"{snp_id} passed PolyMarker design but has no assay "
                    "in the PolyMarker assay table"
                )

            any_valid = any(
                assay_status_by_id[assay_id]["status"] == "PASS"
                for assay_id in assay_ids
            )

            if any_valid:
                status = "PASS"
                failure_reason = ""
            else:
                status = "FAIL"
                failure_reason = (
                    "no_assay_passed_in_silico_validation"
                )

        snp_rows.append({
            "snp_id": snp_id,
            "status": status,
            "failure_reason": failure_reason,
        })

    write_tsv(
        args.assay_status,
        [
            "assay_id",
            "snp_id",
            "status",
            "failure_reason",
        ],
        assay_rows,
    )

    write_tsv(
        args.assay_status_by_genotype,
        [
            "assay_id",
            "genotype",
            "expected_allele",
            "status",
            "failure_reason",
            "expected_amplicons",
            "target_amplicons",
            "unexpected_amplicons",
            "noncanonical_amplicons",
        ],
        genotype_rows,
    )

    write_tsv(
        args.validation_status,
        [
            "snp_id",
            "status",
            "failure_reason",
        ],
        snp_rows,
    )


if __name__ == "__main__":
    main()
