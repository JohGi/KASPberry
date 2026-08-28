from __future__ import annotations

from pathlib import Path

import pytest

from helpers import (
    copy_case_input,
    prepare_test_config,
    run_snakemake,
)

RULE_NAME = "run_mfeprimer_dimer"

TEST_CASES = [
    {
        "case_name": "strong_intra_assay_dimer",
        "target": "results/.work/kasp/in_silico_validation/dimers.tsv",
        "expected_primer_ids": [
            "snp::2::700::assay::01_A_fw",
            "snp::2::700::assay::01_G_fw",
        ],
    },
]


@pytest.mark.parametrize(
    "case",
    TEST_CASES,
    ids=[case["case_name"] for case in TEST_CASES],
)
def test_run_mfeprimer_dimer_cases(
    tmp_path: Path,
    repo_root: Path,
    integration_cases_dir: Path,
    integration_config_dir: Path,
    case: dict[str, object],
) -> None:
    """Run real MFEprimer dimer detection on synthetic primers."""
    case_dir = integration_cases_dir / RULE_NAME / str(case["case_name"])

    copy_case_input(case_dir / "input", tmp_path)

    configfile = prepare_test_config(
        base_config=integration_config_dir / "base_config.yaml",
        workdir=tmp_path,
        override_config=case_dir / "config_override.yaml",
    )

    run_snakemake(
        repo_root=repo_root,
        workdir=tmp_path,
        target=str(case["target"]),
        configfile=configfile,
        allowed_rules=[RULE_NAME],
    )

    output = tmp_path / str(case["target"])
    assert output.exists(), f"Missing MFEprimer dimer output: {output}"
    assert output.stat().st_size > 0, "MFEprimer did not report a dimer"

    text = output.read_text(encoding="utf-8")
    for primer_id in case["expected_primer_ids"]:
        assert str(primer_id) in text, (
            f"Expected dimer primer was not found in report: {primer_id}"
        )
