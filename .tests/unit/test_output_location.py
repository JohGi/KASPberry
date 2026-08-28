from __future__ import annotations

import sys
from pathlib import Path

from snakemake.utils import validate


REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "workflow" / "scripts"))

from run_info import get_snakemake_working_directory


def test_config_schema_accepts_project_without_output_dir() -> None:
    """The scientific configuration only identifies the project."""
    config = {
        "project": {"name": "Test project"},
        "inputs": {"genotypes": "genotypes.tsv"},
        "snps": {},
    }

    validate(
        config,
        REPO_ROOT / "workflow" / "schemas" / "config.schema.yaml",
    )

    assert config["snps"]["min_block_length"] == 301


def test_snakemake_directory_selects_results_parent(
    tmp_path: Path,
) -> None:
    """Snakemake's native directory option determines the results parent."""
    working_directory = get_snakemake_working_directory(
        ["--directory", str(tmp_path)]
    )

    assert working_directory == tmp_path.resolve()
    assert working_directory / "results" == tmp_path / "results"
