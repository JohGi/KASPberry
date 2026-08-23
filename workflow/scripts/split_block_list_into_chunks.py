#!/usr/bin/env python3
"""Split a block list into fixed-size chunk files."""

import argparse
import logging
from pathlib import Path

from attrs import define


logging.basicConfig(
    level=logging.INFO,
    format="[%(levelname)s] %(message)s",
)
LOGGER = logging.getLogger(__name__)


@define(frozen=True)
class ChunkingConfig:
    """Configuration for block list chunking."""

    input_path: Path
    output_dir: Path
    # chunk_map_path: Path
    chunk_size: int


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Split a block list into fixed-size chunk files."
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Input text file containing one block ID per line.",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Output directory where chunk files will be written.",
    )
    # parser.add_argument(
    #     "--chunk-map",
    #     required=True,
    #     help="Output TSV mapping each block ID to its chunk ID.",
    # )
    parser.add_argument(
        "--chunk-size",
        required=True,
        type=int,
        help="Number of block IDs per chunk.",
    )
    return parser.parse_args()


def build_config(args: argparse.Namespace) -> ChunkingConfig:
    """Build a validated configuration object."""
    if args.chunk_size <= 0:
        raise ValueError("--chunk-size must be a strictly positive integer.")

    return ChunkingConfig(
        input_path=Path(args.input),
        output_dir=Path(args.output_dir),
        # chunk_map_path=Path(args.chunk_map),
        chunk_size=args.chunk_size,
    )


def read_block_ids(input_path: Path) -> list[str]:
    """Read block IDs from a text file."""
    if not input_path.is_file():
        raise FileNotFoundError(f"Input block list not found: {input_path}")

    block_ids: list[str] = []
    with input_path.open(encoding="utf-8") as handle:
        for line in handle:
            block_id = line.strip()
            if block_id:
                block_ids.append(block_id)

    if not block_ids:
        raise ValueError(f"Input block list is empty: {input_path}")

    return block_ids


def reset_output_directory(output_dir: Path) -> None:
    """Create an empty output directory."""
    output_dir.mkdir(parents=True, exist_ok=True)
    for chunk_file in output_dir.glob("*.list"):
        chunk_file.unlink()


def write_chunks(
    block_ids: list[str],
    output_dir: Path,
    # chunk_map_path: Path,
    chunk_size: int,
) -> None:
    """Write chunk files and the block-to-chunk mapping."""
    reset_output_directory(output_dir)
    # chunk_map_path.parent.mkdir(parents=True, exist_ok=True)

    # with chunk_map_path.open("w", encoding="utf-8") as chunk_map_handle:
    for start_index in range(0, len(block_ids), chunk_size):
        chunk_index = start_index // chunk_size
        chunk_id = f"chunk_{chunk_index:04d}"
        chunk_path = output_dir / f"{chunk_id}.list"
        chunk_block_ids = block_ids[start_index:start_index + chunk_size]

        LOGGER.info(
            "Writing %s with %d block(s).",
            chunk_path.name,
            len(chunk_block_ids),
        )

        with chunk_path.open("w", encoding="utf-8") as chunk_handle:
            for block_id in chunk_block_ids:
                chunk_handle.write(f"{block_id}\n")
                # chunk_map_handle.write(f"{block_id}\t{chunk_id}\n")


def main() -> None:
    """Run the script."""
    args = parse_args()
    config = build_config(args)
    block_ids = read_block_ids(config.input_path)
    write_chunks(
        block_ids=block_ids,
        output_dir=config.output_dir,
        # chunk_map_path=config.chunk_map_path,
        chunk_size=config.chunk_size,
    )


if __name__ == "__main__":
    main()
