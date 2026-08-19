from __future__ import annotations

from pathlib import Path

import pytest

from helpers import (
    compare_directories,
    copy_case_input,
    prepare_test_config,
    run_snakemake,
)

RULE_NAME = "extract_all_blocks_fasta"

TEST_CASES = [
    {
        "case_name": "case_basic",
        "target": "results/.work/block_fastas/all_blocks.raw.fasta",
    },
]


@pytest.mark.parametrize(
    "case",
    TEST_CASES,
    ids=[case["case_name"] for case in TEST_CASES],
)
def test_extract_all_blocks_fasta_cases(
    tmp_path: Path,
    repo_root: Path,
    integration_cases_dir: Path,
    integration_config_dir: Path,
    case: dict[str, str],
) -> None:
    """Run integration tests for the extract_all_blocks_fasta rule."""
    case_dir = integration_cases_dir / RULE_NAME / case["case_name"]

    copy_case_input(case_dir / "input", tmp_path)

    configfile = prepare_test_config(
        base_config=integration_config_dir / "base_config.yaml",
        workdir=tmp_path,
    )

    run_snakemake(
        repo_root=repo_root,
        workdir=tmp_path,
        target=case["target"],
        configfile=configfile,
    )

    compare_directories(
        expected_dir=case_dir / "expected" / "results" / ".work" / "block_fastas",
        observed_dir=tmp_path / "results" / ".work" / "block_fastas",
    )
