from __future__ import annotations

from pathlib import Path

import pytest

from helpers import (
    compare_directories,
    copy_case_input,
    prepare_test_config,
    run_snakemake,
)

RULE_NAME = "filter_snps_by_groups"

TEST_CASES = [
    {
        "case_name": "1vs1",
        "target": "results/.work/snps/diagnostic_status.tsv",
        "compare_dir": ".work/snps",
    },
    {
        "case_name": "2vs2",
        "target": "results/.work/snps/diagnostic_status.tsv",
        "compare_dir": ".work/snps",
    },
    {
        "case_name": "1vsAll",
        "target": "results/.work/snps/diagnostic_status.tsv",
        "compare_dir": ".work/snps",
    },
]


@pytest.mark.parametrize(
    "case",
    TEST_CASES,
    ids=[case["case_name"] for case in TEST_CASES],
)
def test_filter_snps_by_groups_cases(
    tmp_path: Path,
    repo_root: Path,
    integration_cases_dir: Path,
    integration_config_dir: Path,
    case: dict[str, str],
) -> None:
    """Run integration tests for SNP group filtering."""

    case_dir = integration_cases_dir / RULE_NAME / case["case_name"]

    # Copy input files into tmp working directory
    copy_case_input(case_dir / "input", tmp_path)

    # Prepare config (base + override)
    configfile = prepare_test_config(
        base_config=integration_config_dir / "base_config.yaml",
        workdir=tmp_path,
    )

    # Run Snakemake
    run_snakemake(
        repo_root=repo_root,
        workdir=tmp_path,
        target=case["target"],
        configfile=configfile,
    )

    # Compare outputs
    compare_directories(
        expected_dir=case_dir / "expected" / "results" / case["compare_dir"],
        observed_dir=tmp_path / "results" / case["compare_dir"],
    )
