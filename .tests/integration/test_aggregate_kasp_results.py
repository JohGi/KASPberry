from __future__ import annotations

import filecmp
from pathlib import Path

from helpers import copy_case_input, prepare_test_config, run_snakemake

RULE_NAME = "aggregate_kasp_results"


def test_aggregate_kasp_results_mixed_outcomes(
    tmp_path: Path,
    repo_root: Path,
    integration_cases_dir: Path,
    integration_config_dir: Path,
) -> None:
    """Run the KASP-mode aggregation integration test."""
    case_dir = integration_cases_dir / RULE_NAME / "mixed_kasp_outcomes"

    copy_case_input(case_dir / "input", tmp_path)

    configfile = prepare_test_config(
        base_config=integration_config_dir / "base_config.yaml",
        workdir=tmp_path,
        override_config=case_dir / "config_override.yaml",
    )

    run_snakemake(
        repo_root=repo_root,
        workdir=tmp_path,
        target="results/20_aggregation/kasp_summary.tsv",
        configfile=configfile,
    )

    for relative_path in [
        "results/20_aggregation/kasp_summary.tsv",
        "results/19_validated_assays/validated_snps.vcf",
        "results/assay_summary.tsv",
        "results/validated_assays.tsv",
        "results/primers_to_order.tsv",
    ]:
        assert filecmp.cmp(
            case_dir / "expected" / relative_path,
            tmp_path / relative_path,
            shallow=False,
        ), f"Files differ: {relative_path}"
