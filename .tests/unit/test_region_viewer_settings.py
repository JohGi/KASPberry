from __future__ import annotations

import sys
from pathlib import Path

import yaml


REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "workflow" / "scripts"))

from region_viewer.settings import read_analysis_settings


def write_config(tmp_path: Path, config: dict) -> Path:
    """Write a config YAML and its genotype table for settings tests."""
    genotype_path = tmp_path / "genotypes.tsv"
    genotype_path.write_text(
        "genotype\tgroup\tgenome_fasta\n"
        "alpha\tgroup_a\talpha_genome.fa\n"
        "beta\tgroup_a\t\n"
        "gamma\tgroup_b\tgamma_genome.fa\n",
        encoding="utf-8",
    )
    config["inputs"] = {"genotypes": str(genotype_path)}

    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.safe_dump(config), encoding="utf-8")
    return config_path


def test_read_analysis_settings_for_snps_mode(tmp_path: Path) -> None:
    config_path = write_config(
        tmp_path,
        {
            "snps": {
                "min_block_length": 450,
                "min_snp_flank": 50,
                "repeat_masking": {"library": "references/custom_repeats.fa"},
            },
            "advanced": {
                "sibeliaz": {"extra_options": ["--some-option", "value"]},
                "alignment": {"mafft_options": ["--adjustdirection"]},
                "batching": {"blocks_per_chunk": 2},
            },
        },
    )

    assert read_analysis_settings(config_path, mode="snps") == {
        "minimum_block_length_bp": 450,
        "minimum_snp_flank_bp": 50,
        "snp_groups": {"group_a": ["alpha", "beta"], "group_b": ["gamma"]},
        "repeat_masking": {
            "simple_repeats_and_low_complexity": True,
            "custom_repeat_library": "references/custom_repeats.fa",
        },
        "advanced_options": {
            "sibeliaz": ["--some-option", "value"],
            "mafft": ["--adjustdirection"],
        },
    }


def test_read_analysis_settings_uses_kasp_defaults_and_omits_empty_options(
    tmp_path: Path,
) -> None:
    config_path = write_config(
        tmp_path,
        {
            "snps": {"min_block_length": 450, "min_snp_flank": 50},
            "kasp": {"polymarker_genomes": 2},
            "advanced": {
                "mfeprimer": {"dimer_extra_options": ["--extra-dimer"]},
            },
        },
    )

    settings = read_analysis_settings(config_path, mode="kasp")

    assert settings["repeat_masking"] == {
        "simple_repeats_and_low_complexity": True,
        "custom_repeat_library": None,
    }
    assert settings["kasp_assay_design"] == {
        "polymarker_subgenomes": 2,
        "genotypes": ["alpha", "gamma"],
        "mfeprimer_min_tm": 50,
        "mfeprimer_dimer_max_dg": -3.5,
    }
    assert settings["advanced_options"] == {
        "mfeprimer_dimer": ["--extra-dimer"],
    }
