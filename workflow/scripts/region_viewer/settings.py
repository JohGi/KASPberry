#!/usr/bin/env python3
# Author: Johanna Girodolle

"""Read display-relevant KASPberry analysis settings."""

from __future__ import annotations

import csv
import logging
from pathlib import Path

import yaml


LOGGER = logging.getLogger(__name__)


MFEPRIMER_DEFAULT_MIN_TM = 50
MFEPRIMER_DEFAULT_MAX_DG = -3.5


def read_analysis_settings(
    config_yaml_path: Path | None,
    mode: str | None = None,
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

    inputs = raw.get("inputs") or {}
    snps = raw.get("snps") or {}
    repeat_masking = snps.get("repeat_masking") or {}
    settings: dict[str, object] = {
        "minimum_block_length_bp": snps.get("min_block_length"),
        "minimum_snp_flank_bp": snps.get("min_snp_flank"),
        "snp_groups": read_genotype_groups(inputs.get("genotypes")),
        "repeat_masking": {
            "simple_repeats_and_low_complexity": True,
            "custom_repeat_library": repeat_masking.get("library") or None,
        },
    }

    advanced = raw.get("advanced") or {}
    advanced_options = {
        "sibeliaz": read_option_list(
            (advanced.get("sibeliaz") or {}).get("extra_options")
        ),
        "mafft": read_option_list(
            (advanced.get("alignment") or {}).get("mafft_options")
        ),
    }

    if mode == "kasp":
        kasp = raw.get("kasp") or {}
        mfeprimer = kasp.get("mfeprimer") or {}
        specificity = mfeprimer.get("specificity") or {}
        dimer = mfeprimer.get("dimer") or {}

        settings["kasp_assay_design"] = {
            "polymarker_subgenomes": kasp.get("polymarker_genomes"),
            "genotypes": read_kasp_genotypes(inputs.get("genotypes")),
            "mfeprimer_min_tm": specificity.get(
                "min_tm", MFEPRIMER_DEFAULT_MIN_TM
            ),
            "mfeprimer_dimer_max_dg": dimer.get(
                "max_dg", MFEPRIMER_DEFAULT_MAX_DG
            ),
        }
        advanced_options.update(
            {
                "mfeprimer_specificity": read_option_list(
                    (advanced.get("mfeprimer") or {}).get(
                        "specificity_extra_options"
                    )
                ),
                "mfeprimer_dimer": read_option_list(
                    (advanced.get("mfeprimer") or {}).get(
                        "dimer_extra_options"
                    )
                ),
                "mfeprimer_hairpin": read_option_list(
                    (advanced.get("mfeprimer") or {}).get(
                        "hairpin_extra_options"
                    )
                ),
            }
        )

    non_empty_advanced_options = {
        name: options
        for name, options in advanced_options.items()
        if options
    }
    if non_empty_advanced_options:
        settings["advanced_options"] = non_empty_advanced_options

    return settings


def read_option_list(value: object) -> list[str]:
    """Return a configured command option list without reinterpreting it."""
    if not isinstance(value, list):
        return []

    return [str(option) for option in value]


def read_genotype_groups(
    genotype_path: str | Path | None,
) -> dict[str, list[str]]:
    """Read genotype groups from the genotype table."""
    groups: dict[str, list[str]] = {}

    if not genotype_path:
        return groups

    try:
        with Path(genotype_path).open(
            newline="",
            encoding="utf-8",
        ) as handle:
            reader = csv.DictReader(handle, delimiter="\t")

            for row in reader:
                genotype = (row.get("genotype") or "").strip()
                group = (row.get("group") or "").strip()

                if group:
                    groups.setdefault(group, []).append(genotype)
    except (OSError, csv.Error):
        LOGGER.warning("Could not read genotype TSV: %s", genotype_path)

    return groups


def read_kasp_genotypes(
    genotype_path: str | Path | None,
) -> list[str]:
    """Read genotypes with a whole-genome FASTA available for KASP design/QC."""
    genotypes: list[str] = []

    if not genotype_path:
        return genotypes

    try:
        with Path(genotype_path).open(
            newline="",
            encoding="utf-8",
        ) as handle:
            reader = csv.DictReader(handle, delimiter="\t")

            for row in reader:
                genotype = (row.get("genotype") or "").strip()
                genome_fasta = (row.get("genome_fasta") or "").strip()

                if genotype and genome_fasta:
                    genotypes.append(genotype)

    except (OSError, csv.Error):
        LOGGER.warning("Could not read genotype TSV: %s", genotype_path)

    return genotypes
