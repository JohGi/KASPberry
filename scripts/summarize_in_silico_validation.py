#!/usr/bin/env python3
"""Summarize MFEprimer validation and retain KASP assays that pass all checks."""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path


DIMER_HEADER = re.compile(r"^Dimer \d+:\s+(\S+)\s+x\s+(\S+)$")
DIMER_VALUES = re.compile(r"Score:\s+(\d+),\s+Delta G = (-?\d+(?:\.\d+)?)")
HAIRPIN_HEADER = re.compile(r"^Hairpin \d+:\s+(\S+)$")
HAIRPIN_VALUES = re.compile(
    r"Score:\s+(\d+),\s+Tm = (-?\d+(?:\.\d+)?) °C, "
    r"Delta G = (-?\d+(?:\.\d+)?)"
)
PRIMER_ASSAY = re.compile(r"^(snp::.+::assay::\d+)_(?:comm_rev|[ACGT]_fw)$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--polymarker-summary", required=True, type=Path)
    parser.add_argument("--snp-positions", required=True, type=Path)
    parser.add_argument("--in-silico-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_spec(path: Path) -> list[dict[str, str]]:
    """Read an MFEprimer .spec.tsv, including a header-only no-hit file."""
    with path.open(encoding="utf-8") as handle:
        lines = [line for line in handle if not line.startswith("#1-based coordinate")]
    if not lines:
        return []
    if lines[0].startswith("#"):
        lines[0] = lines[0][1:]
    return list(csv.DictReader(lines, delimiter="\t"))


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def pair_id(hit: dict[str, str]) -> str:
    fp = hit["fpName"].removesuffix("_fp")
    rp = hit["rpName"].removesuffix("_rp")
    if fp != rp:
        raise ValueError(f"Inconsistent primer-pair names in MFEprimer output: {fp} / {rp}")
    return fp


def assay_from_primer(primer_id: str) -> str | None:
    match = PRIMER_ASSAY.match(primer_id)
    return match.group(1) if match else None


def parse_dimers(path: Path) -> list[tuple[str, str, int, float]]:
    records, pending = [], None
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
                    (pending[0], pending[1], int(match.group(1)), float(match.group(2)))
                )
                pending = None
    return records


def parse_hairpins(path: Path) -> list[tuple[str, int, float, float]]:
    records, pending = [], None
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
                    (pending, int(match.group(1)), float(match.group(2)), float(match.group(3)))
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
    args.output_dir.mkdir(parents=True, exist_ok=True)

    snp_status = read_tsv(args.polymarker_summary / "polymarker_snp_status.tsv")
    by_genotype = read_tsv(
        args.polymarker_summary / "polymarker_snp_status_by_genotype.tsv"
    )
    assays = read_tsv(args.polymarker_summary / "polymarker_assays.tsv")
    positions = read_tsv(args.snp_positions)

    snps = {row["marker_id"]: row for row in snp_status}
    expected = {
        (row["marker_id"], row["genotype"]): row["expected_allele"]
        for row in by_genotype
    }
    positions_by_genotype = {
        (f"snp::{row['block_id']}::{row['aln_pos']}", row["sample"]):
        (row["nt"], int(row["pos_in_source_seq"]))
        for row in positions
    }
    genotypes = sorted({row["genotype"] for row in by_genotype})

    canonical = {
        genotype: read_spec(
            args.in_silico_dir / "specificity" / genotype / "canonical.spec.tsv"
        )
        for genotype in genotypes
    }
    noncanonical = {
        genotype: read_spec(
            args.in_silico_dir / "specificity" / genotype / "noncanonical.spec.tsv"
        )
        for genotype in genotypes
    }

    bad_dimers: dict[str, list[tuple[str, str, int, float]]] = defaultdict(list)
    for primer_a, primer_b, score, dg in parse_dimers(args.in_silico_dir / "dimers.tsv"):
        assay_a, assay_b = assay_from_primer(primer_a), assay_from_primer(primer_b)
        if assay_a and assay_a == assay_b:
            bad_dimers[assay_a].append((primer_a, primer_b, score, dg))

    bad_hairpins: dict[str, list[tuple[str, int, float, float]]] = defaultdict(list)
    for primer, score, tm, dg in parse_hairpins(args.in_silico_dir / "hairpins.tsv"):
        assay_id = assay_from_primer(primer)
        if assay_id:
            bad_hairpins[assay_id].append((primer, score, tm, dg))

    assay_rows, genotype_rows, validated = [], [], []

    for assay in assays:
        assay_id, snp_id = assay["assay_id"], assay["marker_id"]
        specificity_results = []

        for genotype in genotypes:
            allele = expected[(snp_id, genotype)]
            position_allele, position = positions_by_genotype[(snp_id, genotype)]
            if allele != position_allele:
                raise ValueError(
                    f"{snp_id}/{genotype}: expected allele {allele} != "
                    f"snp_positions_long.tsv allele {position_allele}"
                )

            if allele == assay["first_allele"]:
                allele_primer, other = assay["first_primer"].upper(), assay["second_allele"]
            elif allele == assay["second_allele"]:
                allele_primer, other = assay["second_primer"].upper(), assay["first_allele"]
            else:
                raise ValueError(f"{snp_id}/{genotype}: allele {allele} absent from {assay_id}")

            expected_pair = f"{assay_id}_{allele}_common"
            unexpected_pair = f"{assay_id}_{other}_common"
            expected_hits = [
                hit for hit in canonical[genotype] if pair_id(hit) == expected_pair
            ]
            unexpected_hits = [
                hit for hit in canonical[genotype] if pair_id(hit) == unexpected_pair
            ]
            noncanonical_hits = [
                hit for hit in noncanonical[genotype]
                if pair_id(hit).startswith(f"{assay_id}_")
            ]
            target_hits = [
                hit for hit in expected_hits
                if target_hit(
                    hit,
                    snps[snp_id]["target_chromosome"],
                    position,
                    allele_primer,
                )
            ]

            passed = (
                len(expected_hits) == 1
                and len(target_hits) == 1
                and not unexpected_hits
                and not noncanonical_hits
            )
            specificity_results.append(passed)
            genotype_rows.append({
                "assay_id": assay_id,
                "genotype": genotype,
                "expected_allele": allele,
                "expected_amplicons": len(expected_hits),
                "target_amplicons": len(target_hits),
                "unexpected_amplicons": len(unexpected_hits),
                "noncanonical_amplicons": len(noncanonical_hits),
                "specificity_pass": str(passed).lower(),
            })

        specificity_pass = all(specificity_results)
        dimer_pass = not bad_dimers[assay_id]
        hairpin_pass = not bad_hairpins[assay_id]
        validation_pass = specificity_pass and dimer_pass and hairpin_pass

        assay_rows.append({
            "assay_id": assay_id,
            "marker_id": snp_id,
            "specificity_pass": str(specificity_pass).lower(),
            "dimer_pass": str(dimer_pass).lower(),
            "hairpin_pass": str(hairpin_pass).lower(),
            "validation_pass": str(validation_pass).lower(),
            "n_genomes_ok": sum(specificity_results),
            "n_genomes_total": len(genotypes),
            "problematic_dimers": len(bad_dimers[assay_id]),
            "problematic_hairpins": len(bad_hairpins[assay_id]),
        })
        if validation_pass:
            validated.append(assay)

    assay_status = {row["assay_id"]: row for row in assay_rows}
    assays_by_snp: dict[str, list[str]] = defaultdict(list)
    for assay in assays:
        assays_by_snp[assay["marker_id"]].append(assay["assay_id"])

    snp_rows = []
    for row in snp_status:
        assay_ids = assays_by_snp[row["marker_id"]]
        n_valid = sum(
            assay_status[assay_id]["validation_pass"] == "true"
            for assay_id in assay_ids
        )
        snp_rows.append({
            "marker_id": row["marker_id"],
            "block_id": row["block_id"],
            "aln_pos": row["aln_pos"],
            "target_chromosome": row["target_chromosome"],
            "polymarker_pass": row["polymarker_pass"],
            "validation_pass": str(n_valid > 0).lower(),
            "n_assays_tested": len(assay_ids),
            "n_assays_validated": n_valid,
        })

    write_tsv(
        args.output_dir / "assay_status.tsv",
        ["assay_id", "marker_id", "specificity_pass", "dimer_pass", "hairpin_pass",
         "validation_pass", "n_genomes_ok", "n_genomes_total",
         "problematic_dimers", "problematic_hairpins"],
        assay_rows,
    )
    write_tsv(
        args.output_dir / "assay_status_by_genotype.tsv",
        ["assay_id", "genotype", "expected_allele", "expected_amplicons",
         "target_amplicons", "unexpected_amplicons", "noncanonical_amplicons",
         "specificity_pass"],
        genotype_rows,
    )
    write_tsv(
        args.output_dir / "snp_status.tsv",
        ["marker_id", "block_id", "aln_pos", "target_chromosome", "polymarker_pass",
         "validation_pass", "n_assays_tested", "n_assays_validated"],
        snp_rows,
    )
    assay_fields = list(assays[0]) if assays else [
        "assay_id", "marker_id", "first_allele", "second_allele", "first_primer",
        "second_primer", "common_primer", "first_primer_with_tail",
        "second_primer_with_tail", "source_genotypes",
    ]
    write_tsv(args.output_dir / "validated_assays.tsv", assay_fields, validated)


if __name__ == "__main__":
    main()
