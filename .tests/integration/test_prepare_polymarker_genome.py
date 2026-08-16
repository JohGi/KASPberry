from __future__ import annotations

from pathlib import Path

import pytest

from helpers import (
    compare_directories,
    copy_case_input,
    copy_shared_resources,
    prepare_test_config,
    run_snakemake,
)

RULE_NAME = "prepare_polymarker_genome"

TEST_CASES = [
    {
        "case_name": "polyploid",
        "target": "results/15_polymarker_inputs/Test/genome.fasta",
    },
    {
        "case_name": "single_genome",
        "target": "results/15_polymarker_inputs/Test/genome.fasta",
    },
]


@pytest.mark.parametrize(
    "case",
    TEST_CASES,
    ids=[case["case_name"] for case in TEST_CASES],
)
def test_prepare_polymarker_genome_cases(
    tmp_path: Path,
    repo_root: Path,
    integration_cases_dir: Path,
    integration_config_dir: Path,
    integration_resources_dir: Path,
    case: dict[str, str],
) -> None:
    """Run integration tests for the mask_block_chunk rule."""
    case_dir = integration_cases_dir / RULE_NAME / case["case_name"]

    copy_case_input(case_dir / "input", tmp_path)
    copy_shared_resources(integration_resources_dir, tmp_path)

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
    )

    compare_directories(
        expected_dir=case_dir / "expected" / "results" / "15_polymarker_inputs",
        observed_dir=tmp_path / "results" / "15_polymarker_inputs"
    )
