from __future__ import annotations

from pathlib import Path

import pytest

from helpers import (
    compare_directories,
    copy_case_input,
    prepare_test_config,
    run_snakemake,
)

RULE_NAME = "filter_sibeliaz_blocks"

TEST_CASES = [
    {
        "case_name": "conserved_4_not_5",
        "target": "results/.work/blocks/collinear_blocks.gff",
    },
    {
        "case_name": "mixed_strand",
        "target": "results/.work/blocks/collinear_blocks.gff",
    },
    {
        "case_name": "too_short",
        "target": "results/.work/blocks/collinear_blocks.gff",
    },
    {
        "case_name": "default_min_length_boundary",
        "target": "results/.work/blocks/collinear_blocks.gff",
    },
]


@pytest.mark.parametrize(
    "case",
    TEST_CASES,
    ids=[case["case_name"] for case in TEST_CASES],
)
def test_filter_sibeliaz_blocks_cases(
    tmp_path: Path,
    repo_root: Path,
    integration_cases_dir: Path,
    integration_config_dir: Path,
    case: dict[str, str],
) -> None:
    """Run integration tests for the filter_sibeliaz_blocks rule."""
    case_dir = integration_cases_dir / RULE_NAME / case["case_name"]

    copy_case_input(case_dir / "input", tmp_path)

    configfile = prepare_test_config(
        base_config=integration_config_dir / "base_config.yaml",
        workdir=tmp_path,
        override_config=case_dir / "config_override.yaml",
    )

    run_snakemake(
        repo_root=repo_root,
        workdir=tmp_path,
        target=case["target"],
        configfile=configfile,
        allowed_rules=[RULE_NAME],
    )

    compare_directories(
        expected_dir=case_dir / "expected" / "results" / ".work" / "blocks",
        observed_dir=tmp_path / "results" / ".work" / "blocks",
    )
