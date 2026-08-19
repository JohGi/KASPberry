from __future__ import annotations

import filecmp
from pathlib import Path

from helpers import copy_case_input, prepare_test_config, run_snakemake

RULE_NAME = "aggregate_snp_results"


def test_aggregate_snp_results_mixed_status(
    tmp_path: Path,
    repo_root: Path,
    integration_cases_dir: Path,
    integration_config_dir: Path,
) -> None:
    """Run the SNP-mode aggregation integration test."""
    case_dir = integration_cases_dir / RULE_NAME / "mixed_status"

    copy_case_input(case_dir / "input", tmp_path)

    configfile = prepare_test_config(
        base_config=integration_config_dir / "base_config.yaml",
        workdir=tmp_path,
    )

    run_snakemake(
        repo_root=repo_root,
        workdir=tmp_path,
        target="results/snps/snp_summary.tsv",
        configfile=configfile,
    )

    for relative_path in [
        "results/snps/snp_summary.tsv",
        "results/snps/diagnostic_snps.vcf",
    ]:
        assert filecmp.cmp(
            case_dir / "expected" / relative_path,
            tmp_path / relative_path,
            shallow=False,
        ), f"Files differ: {relative_path}"
