from __future__ import annotations

from pathlib import Path

import pytest

from helpers import (
    copy_case_input,
    prepare_test_config,
    run_snakemake,
)

RULE_NAME = "run_mfeprimer_specificity"

TEST_CASES = [
    {
        "case_name": "common_common_amplicon",
        "target": (
            "results/.work/kasp/in_silico_validation/"
            "specificity/Test/noncanonical.spec.tsv"
        ),
        "expected_pair_id": (
            "snp::1::100::assay::01_common_common_reverse_self"
        ),
    },
]


@pytest.mark.parametrize(
    "case",
    TEST_CASES,
    ids=[case["case_name"] for case in TEST_CASES],
)
def test_run_mfeprimer_specificity_cases(
    tmp_path: Path,
    repo_root: Path,
    integration_cases_dir: Path,
    integration_config_dir: Path,
    case: dict[str, str],
) -> None:
    """Run real MFEprimer specificity checks on synthetic inputs."""
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
        allowed_rules=[
            RULE_NAME,
            "index_mfeprimer_genome",
        ],
    )

    output = tmp_path / case["target"]
    assert output.exists(), f"Missing MFEprimer specificity output: {output}"
    assert output.stat().st_size > 0, (
        "MFEprimer did not report the expected noncanonical amplicon"
    )

    text = output.read_text(encoding="utf-8")
    assert case["expected_pair_id"] in text, (
        "Expected common x common primer pair was not found in "
        "MFEprimer specificity output"
    )
