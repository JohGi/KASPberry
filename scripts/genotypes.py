"""Shared reader and validation for the public genotypes.tsv interface."""
from __future__ import annotations

import csv
from pathlib import Path

REQUIRED_COLUMNS = {"genotype", "group", "region_fasta", "genome_fasta", "source_seq", "region_start"}

def read_genotypes(path: str | Path) -> list[dict[str, str]]:
    path = Path(path)
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        columns = set(reader.fieldnames or [])
        missing = REQUIRED_COLUMNS - columns
        if missing:
            raise ValueError(f"Missing required columns in {path}: {sorted(missing)}")
        rows = []
        seen = set()
        for line, row in enumerate(reader, 2):
            record = {key: (row.get(key) or "").strip() for key in REQUIRED_COLUMNS}
            name = record["genotype"]
            if not name:
                raise ValueError(f"Empty genotype in {path} at line {line}")
            if name in seen:
                raise ValueError(f"Duplicate genotype {name!r} in {path} at line {line}")
            if not record["region_fasta"]:
                raise ValueError(f"Missing region_fasta for {name!r} in {path} at line {line}")
            if record["genome_fasta"] and (not record["source_seq"] or not record["region_start"]):
                raise ValueError(f"Partial whole-genome metadata for {name!r} in {path} at line {line}")
            if record["region_start"]:
                try:
                    if int(record["region_start"]) < 1: raise ValueError
                except ValueError as exc:
                    raise ValueError(f"Invalid region_start for {name!r} in {path} at line {line}") from exc
            seen.add(name); rows.append(record)
    if not rows: raise ValueError(f"Genotypes table {path} is empty")
    groups = {row["group"] for row in rows if row["group"]}
    if len(groups) == 1: raise ValueError(f"Diagnostic groups in {path} must have exactly two labels, got {sorted(groups)}")
    if len(groups) > 2: raise ValueError(f"Diagnostic groups in {path} must have exactly two labels, got {sorted(groups)}")
    return rows

def read_annotations(path: str | Path, genotype_names: set[str]) -> dict[str, dict[str, str]]:
    with Path(path).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if set(reader.fieldnames or []) < {"genotype", "track", "gff"}:
            raise ValueError(f"Annotations table {path} must contain genotype, track, gff columns")
        result: dict[str, dict[str, str]] = {}
        for row in reader:
            genotype, track, gff = tuple((row.get(k) or "").strip() for k in ("genotype", "track", "gff"))
            if genotype not in genotype_names: raise ValueError(f"Unknown annotation genotype {genotype!r} in {path}")
            result.setdefault(genotype, {})[track] = gff
    return result
