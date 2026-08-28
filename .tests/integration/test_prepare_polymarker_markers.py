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

RULE_NAME = "prepare_polymarker_markers"

TEST_CASES = [
    {
        "case_name": "case_basic",
        "target": "results/.work/kasp/polymarker_inputs/Test/marker_list.csv",
    },
]


@pytest.mark.parametrize(
    "case",
    TEST_CASES,
    ids=[case["case_name"] for case in TEST_CASES],
)
def test_prepare_polymarker_markers_cases(
    tmp_path: Path,
    repo_root: Path,
    integration_cases_dir: Path,
    integration_config_dir: Path,
    integration_resources_dir: Path,
    case: dict[str, str],
) -> None:
    """Run integration tests for the prepare_polymarker_markers rule."""
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
        allowed_rules=[RULE_NAME],
    )

    compare_directories(
        expected_dir=case_dir / "expected" / "results" / ".work" / "kasp" / "polymarker_inputs",
        observed_dir=tmp_path / "results" / ".work" / "kasp" / "polymarker_inputs",
        ignore_suffixes={".tsv", ".fasta"}
    )
