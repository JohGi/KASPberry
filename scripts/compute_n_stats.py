#!/usr/bin/env python3
# Author: Johanna Girodolle

"""Compute paired N-content statistics between masked and unmasked block alignments."""

from __future__ import annotations

import argparse
import logging
import re
from pathlib import Path

LOGGER = logging.getLogger(__name__)

_COORD_SUFFIX_PATTERN = re.compile(r":\d+-\d+$")


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Compute N-content statistics for one block, comparing "
            "masked and unmasked alignments."
        )
    )
    parser.add_argument(
        "--masked-alignment",
        required=True,
        help="Masked block alignment FASTA file (e.g. 3.aln.fasta).",
    )
    parser.add_argument(
        "--unmasked-alignment",
        required=True,
        help="Unmasked block alignment FASTA file (e.g. 3.aln.fasta).",
    )
    return parser.parse_args()


def extract_block_id(fasta_path: Path) -> str:
    """Extract the block ID from an alignment FASTA filename (e.g. '3' from '3.aln.fasta')."""
    return fasta_path.name.split(".", 1)[0]


def extract_sample_name(header: str) -> str:
    """Extract sample name from a FASTA header, removing coordinate suffix (e.g. ':3337-3871')."""
    first_token = header.strip().split()[0]
    return _COORD_SUFFIX_PATTERN.sub("", first_token)


def parse_fasta(fasta_path: Path) -> dict[str, str]:
    """Parse a FASTA file into a sample-name-to-sequence mapping."""
    sequences: dict[str, str] = {}
    current_sample: str | None = None
    current_parts: list[str] = []

    with fasta_path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_sample is not None:
                    sequences[current_sample] = "".join(current_parts)
                current_sample = extract_sample_name(line[1:])
                current_parts = []
                continue
            current_parts.append(line)

    if current_sample is not None:
        sequences[current_sample] = "".join(current_parts)

    return sequences


def count_n(sequence: str) -> int:
    """Count N/n characters in a sequence."""
    return sequence.count("N") + sequence.count("n")


def remove_gaps(sequence: str) -> str:
    """Return the biological sequence without alignment gaps."""
    return sequence.replace("-", "")


def count_repeat_masked_n(unmasked_seq: str, masked_seq: str) -> int:
    """Count bases newly converted to N by repeat masking."""
    if len(unmasked_seq) != len(masked_seq):
        raise ValueError(
            "Masked and unmasked sequences have different ungapped lengths: "
            f"{len(masked_seq)} vs {len(unmasked_seq)}"
        )

    return sum(
        1
        for before, after in zip(unmasked_seq, masked_seq)
        if before.upper() != "N" and after.upper() == "N"
    )


def compute_pct(count: int, length: int) -> float:
    """Compute a percentage; returns 0.0 if length is zero."""
    return 0.0 if length == 0 else (count / length) * 100


def build_rows(masked_path: Path, unmasked_path: Path) -> list[str]:
    """Build TSV rows comparing masked and unmasked alignments for one block.

    Raises ValueError if a sample is missing from the unmasked alignment.
    """
    block_id = extract_block_id(masked_path)
    masked_seqs = parse_fasta(masked_path)
    unmasked_seqs = parse_fasta(unmasked_path)

    rows: list[str] = []

    for sample in sorted(masked_seqs):
        if sample not in unmasked_seqs:
            raise ValueError(
                f"Block {block_id}, sample '{sample}' found in masked alignment "
                f"but missing from unmasked alignment: {unmasked_path}"
            )

        masked_seq = remove_gaps(masked_seqs[sample])
        unmasked_seq = remove_gaps(unmasked_seqs[sample])

        if len(masked_seq) != len(unmasked_seq):
            raise ValueError(
                f"Block {block_id}, sample '{sample}': masked and unmasked "
                f"sequences have different ungapped lengths "
                f"({len(masked_seq)} vs {len(unmasked_seq)})"
            )

        sequence_length_bp = len(unmasked_seq)

        unmasked_n_count = count_n(unmasked_seq)
        masked_n_count = count_n(masked_seq)
        repeat_masked_n_count = count_repeat_masked_n(
            unmasked_seq,
            masked_seq,
        )

        unmasked_n_pct = compute_pct(
            unmasked_n_count,
            sequence_length_bp,
        )
        masked_n_pct = compute_pct(
            masked_n_count,
            sequence_length_bp,
        )
        repeat_masked_n_pct = compute_pct(
            repeat_masked_n_count,
            sequence_length_bp,
        )

        unmasked_length_bp = sequence_length_bp
        masked_length_bp = sequence_length_bp

        rows.append(
            f"{block_id}\t{sample}\t{masked_length_bp}\t{unmasked_length_bp}"
            f"\t{unmasked_n_count}\t{unmasked_n_pct:.2f}"
            f"\t{masked_n_count}\t{masked_n_pct:.2f}"
            f"\t{repeat_masked_n_count}\t{repeat_masked_n_pct:.2f}"
        )

    return rows


def main() -> None:
    """Run the script."""
    logging.basicConfig(level=logging.INFO)
    args = parse_args()

    for row in build_rows(
        masked_path=Path(args.masked_alignment),
        unmasked_path=Path(args.unmasked_alignment),
    ):
        print(row)


if __name__ == "__main__":
    main()
