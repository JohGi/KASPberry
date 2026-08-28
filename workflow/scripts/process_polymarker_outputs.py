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
from ids import parse_snp_id

TAIL_FIRST = "GAAGGTCGGAGTCAACGGATT"
TAIL_SECOND = "GAAGGTGACCAAGTTCATGCT"
SNP_ALLELES_RE = re.compile(r"\[(?P<first>[ACGT])/(?P<second>[ACGT])\]", re.I)
ALLELE_PRIMER_RE = re.compile(
    r"^(?P<marker>.+)(?P<allele>[ACGT])_(?P<kind>1st|2nd)$",
    re.I,
)


@dataclass(frozen=True)
class Candidate:
    snp_id: str
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
    snp_id: str
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
    parser.add_argument("--snp-status", required=True, type=Path)
    parser.add_argument("--snp-status-by-genotype", required=True, type=Path)
    parser.add_argument("--assays", required=True, type=Path)
    parser.add_argument("--mfe-canonical-pairs", required=True, type=Path)
    parser.add_argument("--mfe-noncanonical-pairs", required=True, type=Path)
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

            snp_id, target_chr, sequence = (x.strip() for x in row)

            block_id, aln_pos = parse_snp_id(snp_id)
            allele_match = SNP_ALLELES_RE.search(sequence)

            if allele_match is None:
                raise ValueError(
                    f"{path}:{line_no}: no [A/G]-style SNP in {snp_id}"
                )

            if snp_id in candidates:
                raise ValueError(
                    f"Duplicate marker ID in {path}: {snp_id}"
                )

            candidates[snp_id] = Candidate(
                snp_id=snp_id,
                block_id=block_id,
                aln_pos=aln_pos,
                target_chromosome=target_chr,
                first_allele=allele_match.group("first").upper(),
                second_allele=allele_match.group("second").upper(),
            )

    if not candidates:
        raise ValueError(f"No markers found in {path}")

    return candidates


def empty_proposal() -> dict[str, str]:
    return {
        "first_allele": "",
        "second_allele": "",
        "first": "",
        "second": "",
        "common": "",
    }


def read_primers_to_order(
    path: Path,
    candidates: dict[str, Candidate],
) -> dict[str, dict[str, str]]:
    proposals = {
        snp_id: empty_proposal()
        for snp_id in candidates
    }

    with path.open(encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, start=1):
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            fields = line.split()
            if len(fields) < 2:
                raise ValueError(
                    f"{path}:{line_no}: expected primer ID and sequence"
                )

            primer_id = fields[0]
            sequence = fields[-1].upper()

            if primer_id.endswith("_common"):
                snp_id = primer_id[:-7]
                role = "common"
                allele = ""
            else:
                match = ALLELE_PRIMER_RE.match(primer_id)
                if match is None:
                    raise ValueError(
                        f"Unrecognized PolyMarker primer ID: {primer_id}"
                    )

                snp_id = match.group("marker")
                allele = match.group("allele").upper()
                role = (
                    "first"
                    if match.group("kind").lower() == "1st"
                    else "second"
                )

            if snp_id not in candidates:
                raise ValueError(
                    f"Unknown marker in {path}: {snp_id}"
                )

            proposal = proposals[snp_id]

            if proposal[role] and proposal[role] != sequence:
                raise ValueError(
                    f"Conflicting {role} primers for {snp_id} in {path}"
                )

            proposal[role] = sequence

            if role == "first":
                proposal["first_allele"] = allele
            elif role == "second":
                proposal["second_allele"] = allele

    for snp_id, proposal in proposals.items():
        if not is_complete(proposal):
            continue

        candidate = candidates[snp_id]

        if (
            proposal["first_allele"] != candidate.first_allele
            or proposal["second_allele"] != candidate.second_allele
        ):
            raise ValueError(
                "PolyMarker allele labels do not match marker list for "
                f"{snp_id} in {path}"
            )

    return proposals


def is_complete(proposal: dict[str, str]) -> bool:
    return bool(
        proposal["first"]
        and proposal["second"]
        and proposal["common"]
    )


def strip_tail(
    sequence: str,
    tail: str,
    snp_id: str,
) -> str:
    if not sequence.startswith(tail):
        raise ValueError(
            f"{snp_id}: allele-specific primer lacks expected KASP tail"
        )

    primer = sequence[len(tail):]

    if not primer:
        raise ValueError(
            f"{snp_id}: empty primer after KASP tail removal"
        )

    return primer


def build_unique_assays(
    candidates: dict[str, Candidate],
    genotypes: list[str],
    proposals: dict[str, dict[str, dict[str, str]]],
) -> tuple[list[Assay], dict[str, int]]:
    assays: list[Assay] = []
    unique_counts: dict[str, int] = {}

    for snp_id in candidates:
        marker_proposals = [
            (genotype, proposals[genotype][snp_id])
            for genotype in genotypes
        ]

        if not all(
            is_complete(proposal)
            for _, proposal in marker_proposals
        ):
            unique_counts[snp_id] = 0
            continue

        grouped: dict[
            tuple[str, str, str],
            GroupedAssay,
        ] = {}

        for genotype, proposal in marker_proposals:
            key = (
                proposal["first"],
                proposal["second"],
                proposal["common"],
            )

            if key not in grouped:
                grouped[key] = {
                    "proposal": proposal,
                    "sources": [],
                }

            grouped[key]["sources"].append(genotype)

        unique_counts[snp_id] = len(grouped)

        for index, item in enumerate(grouped.values(), start=1):
            proposal = item["proposal"]

            assays.append(
                {
                    "assay_id": f"{snp_id}::assay::{index:02d}",
                    "snp_id": snp_id,
                    "first_allele": proposal["first_allele"],
                    "second_allele": proposal["second_allele"],
                    "first_primer": strip_tail(
                        proposal["first"],
                        TAIL_FIRST,
                        snp_id,
                    ),
                    "second_primer": strip_tail(
                        proposal["second"],
                        TAIL_SECOND,
                        snp_id,
                    ),
                    "common_primer": proposal["common"],
                    "first_primer_with_tail": proposal["first"],
                    "second_primer_with_tail": proposal["second"],
                    "source_genotypes": sorted(item["sources"]),
                }
            )

    return assays, unique_counts


def write_status_outputs(
    summary_path: Path,
    by_genotype_path: Path,
    candidates: dict[str, Candidate],
    markers_by_genotype: dict[str, dict[str, Candidate]],
    genotypes: list[str],
    proposals: dict[str, dict[str, dict[str, str]]],
    unique_counts: dict[str, int],
) -> None:
    summary_fields = [
        "snp_id",
        "status",
        "failure_reason",
    ]

    by_genotype_fields = [
        "snp_id",
        "genotype",
        "target_chromosome",
        "expected_allele",
        "status",
    ]

    with (
        summary_path.open(
            "w",
            newline="",
            encoding="utf-8",
        ) as summary_handle,
        by_genotype_path.open(
            "w",
            newline="",
            encoding="utf-8",
        ) as by_genotype_handle,
    ):
        summary_writer = csv.DictWriter(
            summary_handle,
            fieldnames=summary_fields,
            delimiter="\t",
            lineterminator="\n",
        )
        by_genotype_writer = csv.DictWriter(
            by_genotype_handle,
            fieldnames=by_genotype_fields,
            delimiter="\t",
            lineterminator="\n",
        )

        summary_writer.writeheader()
        by_genotype_writer.writeheader()

        for snp_id, candidate in candidates.items():
            n_genomes_ok = sum(
                is_complete(proposals[genotype][snp_id])
                for genotype in genotypes
            )

            passed = n_genomes_ok == len(genotypes)
            status = "PASS" if passed else "FAIL"
            failure_reason = "" if passed else "missing_polymarker_assay"

            summary_writer.writerow(
                {
                    "snp_id": snp_id,
                    "status": status,
                    "failure_reason": failure_reason,
                }
            )

            for genotype in genotypes:
                proposal = proposals[genotype][snp_id]
                passed = is_complete(proposal)

                by_genotype_writer.writerow(
                    {
                        "snp_id": snp_id,
                        "genotype": genotype,
                        "target_chromosome": (
                            markers_by_genotype[genotype][snp_id].target_chromosome
                        ),
                        "expected_allele": (
                            markers_by_genotype[genotype][snp_id].first_allele
                        ),
                        "status": "PASS" if passed else "FAIL",
                    }
                )


def write_assay_outputs(
    args: argparse.Namespace,
    assays: list[Assay],
) -> None:
    assay_fields = [
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
    ]

    with (
        args.assays.open(
            "w",
            newline="",
            encoding="utf-8",
        ) as assay_handle,
        args.mfe_canonical_pairs.open(
            "w",
            encoding="utf-8",
        ) as canonical_pairs,
        args.mfe_noncanonical_pairs.open(
            "w",
            encoding="utf-8",
        ) as noncanonical_pairs,
        args.mfe_with_tails.open(
            "w",
            encoding="utf-8",
        ) as fasta,
    ):
        writer = csv.DictWriter(
            assay_handle,
            fieldnames=assay_fields,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()

        for assay in assays:
            assay_row = dict(assay)
            assay_row["source_genotypes"] = ",".join(
                assay["source_genotypes"]
            )
            writer.writerow(assay_row)

            assay_id = assay["assay_id"]
            first_allele = assay["first_allele"]
            second_allele = assay["second_allele"]
            first_primer = assay["first_primer"]
            second_primer = assay["second_primer"]
            common_primer = assay["common_primer"]

            canonical_pairs.write(
                f"{assay_id}_{first_allele}_common\t"
                f"{first_primer}\t{common_primer}\n"
            )
            canonical_pairs.write(
                f"{assay_id}_{second_allele}_common\t"
                f"{second_primer}\t{common_primer}\n"
            )

            noncanonical_pairs.write(
                f"{assay_id}_{first_allele}_{second_allele}"
                f"_forward_forward\t"
                f"{first_primer}\t{second_primer}\n"
            )
            noncanonical_pairs.write(
                f"{assay_id}_{first_allele}_{first_allele}"
                f"_forward_self\t"
                f"{first_primer}\t{first_primer}\n"
            )
            noncanonical_pairs.write(
                f"{assay_id}_{second_allele}_{second_allele}"
                f"_forward_self\t"
                f"{second_primer}\t{second_primer}\n"
            )

            noncanonical_pairs.write(
                f"{assay_id}_common_common_reverse_self\t"
                f"{common_primer}\t{common_primer}\n"
            )

            fasta.write(
                f">{assay_id}_comm_rev\n"
                f"{common_primer}\n"
            )
            fasta.write(
                f">{assay_id}_{first_allele}_fw\n"
                f"{assay['first_primer_with_tail']}\n"
            )
            fasta.write(
                f">{assay_id}_{second_allele}_fw\n"
                f"{assay['second_primer_with_tail']}\n"
            )


def main() -> None:
    args = parse_args()

    marker_paths = index_by_genotype(
        args.marker_lists,
        "marker list",
    )
    pto_paths = index_by_genotype(
        args.primers_to_order,
        "primers_to_order",
    )

    if set(marker_paths) != set(pto_paths):
        raise ValueError(
            "Genotypes differ between marker lists and PolyMarker outputs"
        )

    genotypes = sorted(marker_paths)

    markers_by_genotype = {
        genotype: read_marker_list(marker_paths[genotype])
        for genotype in genotypes
    }

    reference = markers_by_genotype[genotypes[0]]

    for genotype in genotypes[1:]:
        if (
            set(markers_by_genotype[genotype])
            != set(reference)
        ):
            raise ValueError(
                f"marker_list.csv for {genotype} "
                "does not contain the same SNP IDs"
            )

    proposals = {
        genotype: read_primers_to_order(
            pto_paths[genotype],
            markers_by_genotype[genotype],
        )
        for genotype in genotypes
    }

    assays, unique_counts = build_unique_assays(
        reference,
        genotypes,
        proposals,
    )

    write_status_outputs(
        args.snp_status,
        args.snp_status_by_genotype,
        reference,
        markers_by_genotype,
        genotypes,
        proposals,
        unique_counts,
    )

    write_assay_outputs(
        args,
        assays,
    )


if __name__ == "__main__":
    main()
