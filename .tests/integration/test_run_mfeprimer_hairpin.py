from __future__ import annotations

from pathlib import Path

import pytest

from helpers import (
    copy_case_input,
    prepare_test_config,
    run_snakemake,
)

RULE_NAME = "run_mfeprimer_hairpin"

TEST_CASES = [
    {
        "case_name": "strong_hairpin",
        "target": "results/.work/kasp/in_silico_validation/hairpins.tsv",
        "expected_primer_id": "snp::2::800::assay::01_A_fw",
    },
]


@pytest.mark.parametrize(
    "case",
    TEST_CASES,
    ids=[case["case_name"] for case in TEST_CASES],
)
def test_run_mfeprimer_hairpin_cases(
    tmp_path: Path,
    repo_root: Path,
    integration_cases_dir: Path,
    integration_config_dir: Path,
    case: dict[str, str],
) -> None:
    """Run real MFEprimer hairpin detection on a synthetic primer."""
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
        allowed_rules=[RULE_NAME],
    )

    output = tmp_path / case["target"]
    assert output.exists(), f"Missing MFEprimer hairpin output: {output}"
    assert output.stat().st_size > 0, "MFEprimer did not report a hairpin"

    text = output.read_text(encoding="utf-8")
    assert case["expected_primer_id"] in text, (
        "Expected hairpin primer was not found in MFEprimer report"
    )
