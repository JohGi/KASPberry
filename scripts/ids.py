import re


SNP_ID_RE = re.compile(
    r"^snp::(?P<block_id>.+)::(?P<aln_pos>\d+)$"
)


def make_snp_id(block_id: str, aln_pos: str | int) -> str:
    """Build the canonical KASPberry SNP identifier."""
    return f"snp::{block_id}::{aln_pos}"


def parse_snp_id(snp_id: str) -> tuple[str, int]:
    """Parse a canonical KASPberry SNP identifier."""
    match = SNP_ID_RE.fullmatch(snp_id)

    if match is None:
        raise ValueError(f"Invalid SNP ID: {snp_id}")

    return (
        match.group("block_id"),
        int(match.group("aln_pos")),
    )
