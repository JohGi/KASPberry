import importlib.util
from pathlib import Path

import pytest


SCRIPT_PATH = (
    Path(__file__).parents[2]
    / "workflow"
    / "scripts"
    / "summarize_in_silico_validation.py"
)

spec = importlib.util.spec_from_file_location(
    "summarize_in_silico_validation",
    SCRIPT_PATH,
)
module = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(module)

pair_id = module.pair_id


@pytest.mark.parametrize(
    ("fp_name", "rp_name"),
    [
        (
            "snp::12041::343::assay::01_T_common_fp",
            "snp::12041::343::assay::01_T_common_rp",
        ),
        (
            "snp::12041::343::assay::01_T_common_rp",
            "snp::12041::343::assay::01_T_common_fp",
        ),
    ],
)
def test_pair_id_accepts_both_primer_orientations(
    fp_name: str,
    rp_name: str,
) -> None:
    hit = {
        "fpName": fp_name,
        "rpName": rp_name,
    }

    assert pair_id(hit) == "snp::12041::343::assay::01_T_common"


def test_pair_id_rejects_different_pairs() -> None:
    hit = {
        "fpName": "snp::12041::343::assay::01_T_common_fp",
        "rpName": "snp::12041::343::assay::02_T_common_rp",
    }

    with pytest.raises(ValueError, match="Inconsistent primer-pair names"):
        pair_id(hit)
