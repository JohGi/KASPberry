from __future__ import annotations

from pathlib import Path

import pytest

from helpers import (
    compare_directories,
    copy_case_input,
    prepare_test_config,
    run_snakemake,
)

RULE_NAME = "map_snp_coordinates"

TEST_CASES = [
    {
        "case_name": "with_gaps",
        "target": "results/09_snp_positions/snp_positions_long.tsv",
    },
    {
        "case_name": "no_source_seq_offset",
        "target": "results/09_snp_positions/snp_positions_long.tsv",
    },
    {
        "case_name": "vcf_aln_different_order",
        "target": "results/09_snp_positions/snp_positions_long.tsv",
    },
    {
        "case_name": "2snps",
        "target": "results/09_snp_positions/snp_positions_long.tsv",
    },
    {
        "case_name": "2blocks",
        "target": "results/09_snp_positions/snp_positions_long.tsv",
    },
]



@pytest.mark.parametrize(
    "case",
    TEST_CASES,
    ids=[case["case_name"] for case in TEST_CASES],
)
def test_map_snp_coordinates_cases(
    tmp_path: Path,
    repo_root: Path,
    integration_cases_dir: Path,
    integration_config_dir: Path,
    case: dict[str, str],
) -> None:
    """Run integration tests for the map_snp_coordinates rule."""
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
        expected_dir=case_dir / "expected" / "results" / "09_snp_positions",
        observed_dir=tmp_path / "results" / "09_snp_positions",
    )
