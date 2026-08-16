#!/usr/bin/env python3
"""Consolidate PolyMarker outputs across genotypes and prepare MFEprimer inputs.

For every candidate SNP, this script records whether PolyMarker produced a complete
KASP assay in each genotype. SNPs designable in all genotypes are retained, all
proposed assays are gathered and exact sequence duplicates are collapsed, then
MFEprimer input files are written for the unique assays.
"""

from __future__ import annotations

import argparse
import csv
import re
from dataclasses import dataclass
from pathlib import Path
from typing import TypedDict


TAIL_FIRST = "GAAGGTCGGAGTCAACGGATT"
TAIL_SECOND = "GAAGGTGACCAAGTTCATGCT"
SNP_ID_RE = re.compile(r"^snp::(?P<block>.+)::(?P<pos>\d+)$")
SNP_ALLELES_RE = re.compile(r"\[(?P<first>[ACGT])/(?P<second>[ACGT])\]", re.I)
ALLELE_PRIMER_RE = re.compile(r"^(?P<marker>.+)(?P<allele>[ACGT])_(?P<kind>1st|2nd)$", re.I)
DETAIL_FIELDS = ("chromosome", "SNP_type", "primer_type", "orientation", "total_contigs", "errors", "repetitive", "total_hits")


@dataclass(frozen=True)
class Candidate:
    marker_id: str
    block_id: str
    aln_pos: int
    target_chromosome: str
    first_allele: str
    second_allele: str


class GroupedAssay(TypedDict):
    proposal: dict[str, str]
    sources: list[str]


class Assay(TypedDict):
    assay_id: str
    marker_id: str
    block_id: str
    aln_pos: int
    first_allele: str
    second_allele: str
    first_primer: str
    second_primer: str
    common_primer: str
    first_primer_with_tail: str
    second_primer_with_tail: str
    source_genotypes: list[str]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--marker-lists", nargs="+", required=True, type=Path)
    parser.add_argument("--primers-to-order", nargs="+", required=True, type=Path)
    parser.add_argument("--primers-csv", nargs="+", required=True, type=Path)
    parser.add_argument("--snp-status", required=True, type=Path)
    parser.add_argument("--snp-status-long", required=True, type=Path)
    parser.add_argument("--assays", required=True, type=Path)
    parser.add_argument("--mfe-pairs", required=True, type=Path)
    parser.add_argument("--mfe-control-pairs", required=True, type=Path)
    parser.add_argument("--mfe-with-tails", required=True, type=Path)
    return parser.parse_args()


def index_by_genotype(paths: list[Path], label: str) -> dict[str, Path]:
    result: dict[str, Path] = {}
    for path in paths:
        genotype = path.parent.name
        if genotype in result:
            raise ValueError(f"Duplicate {label} for genotype {genotype}")
        result[genotype] = path
    return result


def read_marker_list(path: Path) -> dict[str, Candidate]:
    candidates: dict[str, Candidate] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for line_no, row in enumerate(csv.reader(handle), start=1):
            if not row:
                continue
            if len(row) != 3:
                raise ValueError(f"{path}:{line_no}: expected 3 columns")

            marker_id, target_chr, sequence = (x.strip() for x in row)
            id_match = SNP_ID_RE.match(marker_id)
            allele_match = SNP_ALLELES_RE.search(sequence)
            if id_match is None:
                raise ValueError(f"Unexpected marker ID: {marker_id}")
            if allele_match is None:
                raise ValueError(f"{path}:{line_no}: no [A/G]-style SNP in {marker_id}")
            if marker_id in candidates:
                raise ValueError(f"Duplicate marker ID in {path}: {marker_id}")

            candidates[marker_id] = Candidate(
                marker_id=marker_id,
                block_id=id_match.group("block"),
                aln_pos=int(id_match.group("pos")),
                target_chromosome=target_chr,
                first_allele=allele_match.group("first").upper(),
                second_allele=allele_match.group("second").upper(),
            )

    if not candidates:
        raise ValueError(f"No markers found in {path}")
    return candidates


def empty_proposal() -> dict[str, str]:
    return {"first_allele": "", "second_allele": "", "first": "", "second": "", "common": ""}


def read_primers_to_order(path: Path, candidates: dict[str, Candidate]) -> dict[str, dict[str, str]]:
    proposals = {marker_id: empty_proposal() for marker_id in candidates}

    with path.open(encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 2:
                raise ValueError(f"{path}:{line_no}: expected primer ID and sequence")

            primer_id, sequence = fields[0], fields[-1].upper()
            if primer_id.endswith("_common"):
                marker_id = primer_id[:-7]
                role = "common"
                allele = ""
            else:
                match = ALLELE_PRIMER_RE.match(primer_id)
                if match is None:
                    raise ValueError(f"Unrecognized PolyMarker primer ID: {primer_id}")
                marker_id = match.group("marker")
                allele = match.group("allele").upper()
                role = "first" if match.group("kind").lower() == "1st" else "second"

            if marker_id not in candidates:
                raise ValueError(f"Unknown marker in {path}: {marker_id}")
            proposal = proposals[marker_id]
            if proposal[role] and proposal[role] != sequence:
                raise ValueError(f"Conflicting {role} primers for {marker_id} in {path}")
            proposal[role] = sequence
            if role == "first":
                proposal["first_allele"] = allele
            elif role == "second":
                proposal["second_allele"] = allele

    for marker_id, proposal in proposals.items():
        if not is_complete(proposal):
            continue
        candidate = candidates[marker_id]
        if proposal["first_allele"] != candidate.first_allele or proposal["second_allele"] != candidate.second_allele:
            raise ValueError(f"PolyMarker allele labels do not match marker list for {marker_id} in {path}")

    return proposals


def is_complete(proposal: dict[str, str]) -> bool:
    return bool(proposal["first"] and proposal["second"] and proposal["common"])


def missing_roles(proposal: dict[str, str]) -> list[str]:
    return [role for role in ("first", "second", "common") if not proposal[role]]


def strip_tail(sequence: str, tail: str, marker_id: str) -> str:
    if not sequence.startswith(tail):
        raise ValueError(f"{marker_id}: allele-specific primer lacks expected KASP tail")
    primer = sequence[len(tail):]
    if not primer:
        raise ValueError(f"{marker_id}: empty primer after KASP tail removal")
    return primer


def read_primer_details(path: Path, candidates: dict[str, Candidate]) -> dict[str, dict[str, str]]:
    details: dict[str, dict[str, str]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or "Marker" not in reader.fieldnames:
            return details
        for row in reader:
            marker_id = (row.get("Marker") or "").strip()
            if marker_id in candidates:
                details[marker_id] = {field: (row.get(field) or "").strip() for field in DETAIL_FIELDS}
    return details


def build_unique_assays(
    candidates: dict[str, Candidate],
    genotypes: list[str],
    proposals: dict[str, dict[str, dict[str, str]]],
) -> tuple[list[Assay], dict[str, int]]:
    assays: list[Assay] = []
    unique_counts: dict[str, int] = {}

    for marker_id, candidate in candidates.items():
        marker_proposals = [(g, proposals[g][marker_id]) for g in genotypes]
        if not all(is_complete(p) for _, p in marker_proposals):
            unique_counts[marker_id] = 0
            continue

        grouped: dict[tuple[str, str, str], GroupedAssay] = {}
        for genotype, proposal in marker_proposals:
            key = (proposal["first"], proposal["second"], proposal["common"])
            if key not in grouped:
                grouped[key] = {"proposal": proposal, "sources": []}
            grouped[key]["sources"].append(genotype)

        unique_counts[marker_id] = len(grouped)
        for index, item in enumerate(grouped.values(), start=1):
            proposal = item["proposal"]
            assays.append({
                "assay_id": f"{marker_id}::assay::{index:02d}",
                "marker_id": marker_id,
                "block_id": candidate.block_id,
                "aln_pos": candidate.aln_pos,
                "first_allele": proposal["first_allele"],
                "second_allele": proposal["second_allele"],
                "first_primer": strip_tail(proposal["first"], TAIL_FIRST, marker_id),
                "second_primer": strip_tail(proposal["second"], TAIL_SECOND, marker_id),
                "common_primer": proposal["common"],
                "first_primer_with_tail": proposal["first"],
                "second_primer_with_tail": proposal["second"],
                "source_genotypes": sorted(item["sources"]),
            })

    return assays, unique_counts


def write_status_outputs(
    summary_path: Path,
    long_path: Path,
    candidates: dict[str, Candidate],
    genotypes: list[str],
    proposals: dict[str, dict[str, dict[str, str]]],
    details: dict[str, dict[str, dict[str, str]]],
    unique_counts: dict[str, int],
) -> None:
    summary_fields = ["marker_id", "block_id", "aln_pos", "target_chromosome", "n_genomes_ok", "n_genomes_total", "polymarker_pass", "failed_genotypes", "n_assays_proposed", "n_unique_assays_to_test"]
    long_fields = ["marker_id", "block_id", "aln_pos", "target_chromosome", "genotype", "polymarker_pass", "failure_reason", *DETAIL_FIELDS]

    with summary_path.open("w", newline="", encoding="utf-8") as summary_handle, long_path.open("w", newline="", encoding="utf-8") as long_handle:
        summary = csv.DictWriter(summary_handle, fieldnames=summary_fields, delimiter="\t", lineterminator="\n")
        long = csv.DictWriter(long_handle, fieldnames=long_fields, delimiter="\t", lineterminator="\n")
        summary.writeheader()
        long.writeheader()

        for marker_id, candidate in candidates.items():
            ok = [g for g in genotypes if is_complete(proposals[g][marker_id])]
            failed = [g for g in genotypes if g not in ok]
            summary.writerow({
                "marker_id": marker_id,
                "block_id": candidate.block_id,
                "aln_pos": candidate.aln_pos,
                "target_chromosome": candidate.target_chromosome,
                "n_genomes_ok": len(ok),
                "n_genomes_total": len(genotypes),
                "polymarker_pass": str(not failed).lower(),
                "failed_genotypes": ",".join(failed),
                "n_assays_proposed": len(ok),
                "n_unique_assays_to_test": unique_counts[marker_id],
            })

            for genotype in genotypes:
                proposal = proposals[genotype][marker_id]
                detail = details[genotype].get(marker_id, {})
                missing = missing_roles(proposal)
                long.writerow({
                    "marker_id": marker_id,
                    "block_id": candidate.block_id,
                    "aln_pos": candidate.aln_pos,
                    "target_chromosome": candidate.target_chromosome,
                    "genotype": genotype,
                    "polymarker_pass": str(not missing).lower(),
                    "failure_reason": "" if not missing else "missing:" + ",".join(missing),
                    **{field: detail.get(field, "") for field in DETAIL_FIELDS},
                })


def write_assay_outputs(args: argparse.Namespace, assays: list[Assay]) -> None:
    assay_fields = ["assay_id", "marker_id", "block_id", "aln_pos", "first_allele", "second_allele", "first_primer", "second_primer", "common_primer", "first_primer_with_tail", "second_primer_with_tail", "source_genotypes"]

    with args.assays.open("w", newline="", encoding="utf-8") as assay_handle, args.mfe_pairs.open("w", encoding="utf-8") as pairs, args.mfe_control_pairs.open("w", encoding="utf-8") as controls, args.mfe_with_tails.open("w", encoding="utf-8") as fasta:
        writer = csv.DictWriter(assay_handle, fieldnames=assay_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()

        for assay in assays:
            assay_row = dict(assay)
            assay_row["source_genotypes"] = ",".join(assay["source_genotypes"])
            writer.writerow(assay_row)

            aid = assay["assay_id"]
            a1, a2 = assay["first_allele"], assay["second_allele"]
            p1, p2, common = assay["first_primer"], assay["second_primer"], assay["common_primer"]
            pairs.write(f"{aid}_{a1}_common\t{p1}\t{common}\n")
            pairs.write(f"{aid}_{a2}_common\t{p2}\t{common}\n")
            controls.write(f"{aid}_{a1}_{a2}_forward_forward\t{p1}\t{p2}\n")
            controls.write(f"{aid}_{a1}_{a1}_forward_self\t{p1}\t{p1}\n")
            controls.write(f"{aid}_{a2}_{a2}_forward_self\t{p2}\t{p2}\n")
            fasta.write(f">{aid}_comm_rev\n{common}\n")
            fasta.write(f">{aid}_{a1}_fw\n{assay['first_primer_with_tail']}\n")
            fasta.write(f">{aid}_{a2}_fw\n{assay['second_primer_with_tail']}\n")


def main() -> None:
    args = parse_args()

    marker_paths = index_by_genotype(args.marker_lists, "marker list")
    pto_paths = index_by_genotype(args.primers_to_order, "primers_to_order")
    primer_paths = index_by_genotype(args.primers_csv, "primers.csv")
    if not (set(marker_paths) == set(pto_paths) == set(primer_paths)):
        raise ValueError("Genotypes differ between marker lists and PolyMarker outputs")

    genotypes = sorted(marker_paths)
    markers_by_genotype = {g: read_marker_list(marker_paths[g]) for g in genotypes}
    reference = markers_by_genotype[genotypes[0]]
    for genotype in genotypes[1:]:
        if set(markers_by_genotype[genotype]) != set(reference):
            raise ValueError(f"marker_list.csv for {genotype} does not contain the same SNP IDs")

    proposals = {g: read_primers_to_order(pto_paths[g], markers_by_genotype[g]) for g in genotypes}
    details = {g: read_primer_details(primer_paths[g], markers_by_genotype[g]) for g in genotypes}
    assays, unique_counts = build_unique_assays(reference, genotypes, proposals)

    write_status_outputs(args.snp_status, args.snp_status_long, reference, genotypes, proposals, details, unique_counts)
    write_assay_outputs(args, assays)


if __name__ == "__main__":
    main()
