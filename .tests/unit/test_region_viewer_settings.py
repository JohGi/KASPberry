from __future__ import annotations

import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "workflow" / "scripts"))

from region_viewer.settings import read_analysis_settings


def write_analysis_settings(tmp_path: Path, settings: dict) -> Path:
    """Write the resolved settings JSON consumed by the viewer."""
    settings_path = tmp_path / "analysis_settings.json"
    settings_path.write_text(
        json.dumps(settings),
        encoding="utf-8",
    )
    return settings_path


def test_read_analysis_settings_for_snps_mode(tmp_path: Path) -> None:
    settings_path = write_analysis_settings(
        tmp_path,
        {
            "minimum_block_length_bp": 450,
            "snp_groups": {
                "group_a": ["alpha", "beta"],
                "group_b": ["gamma"],
            },
            "repeat_masking": {
                "simple_repeats_and_low_complexity": True,
                "custom_repeat_library": "references/custom_repeats.fa",
            },
            "advanced_options": {
                "sibeliaz": ["--some-option", "value"],
                "mafft": ["--adjustdirection"],
            },
        },
    )

    assert read_analysis_settings(settings_path) == {
        "minimum_block_length_bp": 450,
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


def test_read_analysis_settings_reads_schema_resolved_kasp_defaults(
    tmp_path: Path,
) -> None:
    settings_path = write_analysis_settings(
        tmp_path,
        {
            "minimum_block_length_bp": 301,
            "snp_groups": {},
            "repeat_masking": {
                "simple_repeats_and_low_complexity": True,
                "custom_repeat_library": None,
            },
            "kasp_assay_design": {
                "polymarker_subgenomes": 2,
                "genotypes": ["alpha", "gamma"],
                "mfeprimer_min_tm": 50,
                "mfeprimer_dimer_max_dg": -3.5,
            },
            "advanced_options": {
                "mfeprimer_dimer": ["--extra-dimer"],
            },
        },
    )

    settings = read_analysis_settings(settings_path)

    assert settings["minimum_block_length_bp"] == 301
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
