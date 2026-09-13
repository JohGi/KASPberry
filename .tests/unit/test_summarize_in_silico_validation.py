import importlib.util
from pathlib import Path

import pytest
import sys


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
sys.modules[spec.name] = module
assert spec.loader is not None
spec.loader.exec_module(module)

pair_id = module.pair_id
target_hit = module.target_hit
ValidationContext = module.ValidationContext
evaluate_assay_in_genotype = module.evaluate_assay_in_genotype

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
        target_alias="2A__001",
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
        target_alias="2A__001",
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
        target_alias="2A__001",
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
            target_alias="2A__001",
            target_position=343,
        )


def test_target_hit_rejects_different_sequence_with_same_polymarker_prefix() -> None:
    hit = {
        "chrom": "2A__002",
        "fpName": "snp::12041::343::assay::01_T_common_fp",
        "rpName": "snp::12041::343::assay::01_T_common_rp",
        "fpEnd": "343",
        "rpStart": "400",
    }

    assert not target_hit(
        hit,
        target_alias="2A__001",
        target_position=343,
    )


@pytest.mark.parametrize(
    ("pair_suffix", "count_column"),
    [
        (
            "A_A_allele_self",
            "allele1/allele1_amplicons",
        ),
        (
            "G_G_allele_self",
            "allele2/allele2_amplicons",
        ),
        (
            "common_self",
            "common/common_amplicons",
        ),
    ],
)
def test_noncanonical_pair_counts_are_reported_separately(
    pair_suffix: str,
    count_column: str,
) -> None:
    assay_id = "snp::12041::343::assay::01"
    snp_id = "snp::12041::343"

    assay = {
        "assay_id": assay_id,
        "snp_id": snp_id,
        "first_allele": "A",
        "second_allele": "G",
    }

    target_hit_row = {
        "chrom": "2A__001",
        "fpName": f"{assay_id}_A_common_fp",
        "rpName": f"{assay_id}_A_common_rp",
        "fpEnd": "343",
        "rpStart": "400",
    }

    noncanonical_hit = {
        "fpName": f"{assay_id}_{pair_suffix}_fp",
        "rpName": f"{assay_id}_{pair_suffix}_rp",
    }

    context = ValidationContext(
        genotypes=["A"],
        design_by_snp={},
        design_by_snp_genotype={
            (snp_id, "A"): {
                "expected_allele": "A",
            },
        },
        positions_by_genotype={
            (snp_id, "A"): ("A", 343),
        },
        source_aliases={
            "A": "2A__001",
        },
        canonical={
            "A": [target_hit_row],
        },
        noncanonical={
            "A": [noncanonical_hit],
        },
    )

    row, failure_reasons = evaluate_assay_in_genotype(
        assay,
        "A",
        context,
    )

    assert row["target_amplicons"] == 1
    assert row["allele1/common_amplicons"] == 1
    assert row["allele2/common_amplicons"] == 0
    assert row[count_column] == 1

    other_noncanonical_columns = {
        "allele1/allele2_amplicons",
        "allele1/allele1_amplicons",
        "allele2/allele2_amplicons",
        "common/common_amplicons",
    } - {count_column}

    for column in other_noncanonical_columns:
        assert row[column] == 0

    assert row["status"] == "FAIL"
    assert row["failure_reason"] == "unexpected_amplicon"
    assert failure_reasons == ["unexpected_amplicon"]
