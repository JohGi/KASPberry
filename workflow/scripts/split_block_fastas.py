#!/usr/bin/env python3
"""Split a global block FASTA file into one FASTA per block."""

import argparse
from pathlib import Path
from typing import Iterator


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Split a global block FASTA into one FASTA per block."
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Input FASTA containing all extracted block sequences.",
    )
    parser.add_argument(
        "--block-list",
        required=True,
        help="Text file listing expected block IDs, one per line.",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory for per-block FASTA files.",
    )
    return parser.parse_args()


def read_block_ids(path: Path) -> set[str]:
    """Read expected block IDs."""
    with path.open(encoding="utf-8") as handle:
        return {line.strip() for line in handle if line.strip()}


def iter_fasta_records(path: Path) -> Iterator[tuple[str, str]]:
    """Yield FASTA records as (header, sequence)."""
    header = None
    sequence_lines: list[str] = []

    with path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence_lines)
                header = line[1:]
                sequence_lines = []
            else:
                sequence_lines.append(line)

    if header is not None:
        yield header, "".join(sequence_lines)


def split_block_prefixed_header(header: str) -> tuple[str, str]:
    """Split a FASTA header into block ID and cleaned sequence header.

    Expected format:
        <block_id>__<sequence_header>

    Example:
        30090__Belalur:20506376-20515003
    returns:
        ("30090", "Belalur:20506376-20515003")
    """
    if "__" not in header:
        raise ValueError(
            f"Invalid FASTA header {header!r}: expected '<block_id>__<sequence_header>'."
        )

    block_id, sequence_header = header.split("__", 1)

    if not block_id:
        raise ValueError(f"Invalid FASTA header {header!r}: empty block ID.")

    if not sequence_header:
        raise ValueError(f"Invalid FASTA header {header!r}: empty sequence header.")

    return block_id, sequence_header


def write_block_fastas(
    input_fasta: Path,
    expected_block_ids: set[str],
    outdir: Path,
) -> None:
    """Write one FASTA file per block with cleaned sequence headers."""
    outdir.mkdir(parents=True, exist_ok=True)

    for block_id in sorted(expected_block_ids):
        output_path = outdir / f"{block_id}.fasta"
        if output_path.exists():
            output_path.unlink()

    seen_block_ids: set[str] = set()

    for header, sequence in iter_fasta_records(input_fasta):
        block_id, cleaned_header = split_block_prefixed_header(header)

        if block_id not in expected_block_ids:
            raise ValueError(
                f"Unexpected block ID {block_id!r} found in FASTA header {header!r}."
            )

        output_path = outdir / f"{block_id}.fasta"
        with output_path.open("a", encoding="utf-8") as handle:
            handle.write(f">{cleaned_header}\n{sequence}\n")

        seen_block_ids.add(block_id)

    missing_block_ids = expected_block_ids - seen_block_ids
    if missing_block_ids:
        raise ValueError(
            f"No FASTA entries were written for block IDs: {sorted(missing_block_ids)}"
        )


def main() -> None:
    """Run the script."""
    args = parse_args()
    input_fasta = Path(args.input)
    block_list = Path(args.block_list)
    outdir = Path(args.outdir)

    expected_block_ids = read_block_ids(block_list)
    write_block_fastas(
        input_fasta=input_fasta,
        expected_block_ids=expected_block_ids,
        outdir=outdir,
    )


if __name__ == "__main__":
    main()
