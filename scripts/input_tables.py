"""Shared reader and validation for the public genotypes.tsv interface."""
from __future__ import annotations

import csv
from pathlib import Path

GENOTYPE_COLUMNS = (
    "genotype",
    "group",
    "region_fasta",
    "genome_fasta",
    "source_seq",
    "region_start",
)

def read_genotypes(path: str | Path) -> list[dict[str, str]]:
    """Read and normalize genotype records from a TSV file."""
    path = Path(path)
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [
            {key: (row.get(key) or "").strip() for key in GENOTYPE_COLUMNS}
            for row in reader
        ]

def read_annotations(path: str | Path) -> dict[str, dict[str, str]]:
    """Read annotations into the nested mapping consumed by the viewer."""
    with Path(path).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        result: dict[str, dict[str, str]] = {}
        for row in reader:
            genotype, track, gff = tuple((row.get(k) or "").strip() for k in ("genotype", "track", "gff"))
            result.setdefault(genotype, {})[track] = gff
    return result
