#!/usr/bin/env python3
# Author: Johanna Girodolle

"""Read display-relevant KASPberry analysis settings."""

from __future__ import annotations

import csv
import logging
from pathlib import Path

import yaml


LOGGER = logging.getLogger(__name__)


def read_analysis_settings(
    config_yaml_path: Path | None,
) -> dict[str, object]:
    """Extract display-relevant analysis settings from the pipeline config."""
    if config_yaml_path is None:
        return {}

    try:
        raw = yaml.safe_load(
            config_yaml_path.read_text(encoding="utf-8")
        ) or {}
    except (OSError, yaml.YAMLError):
        LOGGER.warning("Could not read config YAML: %s", config_yaml_path)
        return {}

    snps = raw.get("snps") or {}

    return {
        "minimum_block_length_bp": snps.get("min_block_length"),
        "minimum_snp_flank_bp": snps.get("min_snp_flank"),
        "snp_groups": read_genotype_groups(
            Path(raw["inputs"]["genotypes"])
        ),
    }


def read_genotype_groups(
    genotype_path: Path,
) -> dict[str, list[str]]:
    """Read genotype groups from the genotype table."""
    groups: dict[str, list[str]] = {}

    with genotype_path.open(
        newline="",
        encoding="utf-8",
    ) as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        for row in reader:
            genotype = row["genotype"].strip()
            group = row["group"].strip()

            if group:
                groups.setdefault(group, []).append(genotype)

    return groups
