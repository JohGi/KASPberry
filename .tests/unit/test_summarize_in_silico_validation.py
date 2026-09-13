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
target_hit = module.target_hit

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


def test_target_hit_accepts_allele_primer_on_forward_side() -> None:
    hit = {
        "chrom": "2A__001",
        "fpName": "snp::12041::343::assay::01_T_common_fp",
        "rpName": "snp::12041::343::assay::01_T_common_rp",
        "fpEnd": "343",
        "rpStart": "400",
    }

    assert target_hit(
        hit,
        target_chromosome="2A",
        target_position=343,
    )


def test_target_hit_accepts_allele_primer_on_reverse_side() -> None:
    hit = {
        "chrom": "2A__001",
        "fpName": "snp::12041::343::assay::01_T_common_rp",
        "rpName": "snp::12041::343::assay::01_T_common_fp",
        "fpEnd": "400",
        "rpStart": "343",
    }

    assert target_hit(
        hit,
        target_chromosome="2A",
        target_position=343,
    )


def test_target_hit_rejects_wrong_reverse_three_prime_position() -> None:
    hit = {
        "chrom": "2A__001",
        "fpName": "snp::12041::343::assay::01_T_common_rp",
        "rpName": "snp::12041::343::assay::01_T_common_fp",
        "fpEnd": "400",
        "rpStart": "344",
    }

    assert not target_hit(
        hit,
        target_chromosome="2A",
        target_position=343,
    )


def test_target_hit_rejects_invalid_pair_members() -> None:
    hit = {
        "chrom": "2A__001",
        "fpName": "snp::12041::343::assay::01_T_common_fp",
        "rpName": "snp::12041::343::assay::01_T_common_fp",
        "fpEnd": "343",
        "rpStart": "343",
    }

    with pytest.raises(
        ValueError,
        match="Expected exactly one _fp and one _rp member",
    ):
        target_hit(
            hit,
            target_chromosome="2A",
            target_position=343,
        )
