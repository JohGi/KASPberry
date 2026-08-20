from __future__ import annotations

from pathlib import Path

from helpers import (
    copy_case_input,
    copy_shared_resources,
    prepare_test_config,
    run_snakemake,
)

RULE_NAME = "compute_masked_block_n_stats"
CASE_NAME = "case_basic"
TARGET = "results/.work/block_stats/masked_block_n_stats.tsv"


def test_compute_masked_block_n_stats(
    tmp_path: Path,
    repo_root: Path,
    integration_cases_dir: Path,
    integration_config_dir: Path,
    integration_resources_dir: Path,
) -> None:
    """Check ungapped N-content and newly repeat-masked base percentages."""
    case_dir = integration_cases_dir / RULE_NAME / CASE_NAME

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
        target=TARGET,
        configfile=configfile,
    )

    expected = case_dir / "expected" / TARGET
    observed = tmp_path / TARGET

    assert observed.read_text(encoding="utf-8") == expected.read_text(
        encoding="utf-8"
    )
