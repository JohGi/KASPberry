#!/usr/bin/env python3
"""Command-line interface for KASPberry."""

from __future__ import annotations
from snakemake_interface_common.exceptions import WorkflowError

import argparse
import subprocess
import sys
from pathlib import Path

import yaml

from validate_inputs import validate_inputs


REPO_ROOT = Path(__file__).resolve().parents[1]
SNAKEFILE = REPO_ROOT / "Snakefile"
SCHEMA_DIR = REPO_ROOT / "workflow" / "schemas"


def parse_args() -> tuple[argparse.Namespace, list[str]]:
    """Parse KASPberry arguments and preserve extra arguments for Snakemake."""
    parser = argparse.ArgumentParser(
        prog="kaspberry",
        description="Discover diagnostic SNPs and prioritize KASP assay candidates.",
    )

    subparsers = parser.add_subparsers(
        dest="mode",
        required=True,
    )

    for mode in ("snps", "kasp"):
        subparser = subparsers.add_parser(
            mode,
            help=f"Run the KASPberry {mode} workflow.",
        )
        subparser.add_argument(
            "-c",
            "--config",
            required=True,
            type=Path,
            help="Path to the KASPberry YAML configuration file.",
        )

    # Arguments unknown to KASPberry are forwarded to Snakemake.
    return parser.parse_known_args()


def load_config(config_path: Path) -> dict:
    """Load a KASPberry YAML configuration file."""
    if not config_path.is_file():
        raise ValueError(f"Configuration file not found: {config_path}")

    with config_path.open("r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle)

    if not isinstance(config, dict):
        raise ValueError(
            f"Configuration file must contain a YAML mapping: {config_path}"
        )

    return config


def run_snakemake(
    mode: str,
    config_path: Path,
    extra_args: list[str],
) -> int:
    """Run the requested Snakemake target."""
    command = [
        "snakemake",
        "--snakefile",
        str(SNAKEFILE),
        "--configfile",
        str(config_path),
        *extra_args,
        mode,
    ]

    return subprocess.run(command).returncode


def main() -> int:
    """Run the KASPberry CLI."""
    args, snakemake_args = parse_args()

    config_path = args.config.resolve()

    try:
        config = load_config(config_path)

        validate_inputs(
            config=config,
            mode=args.mode,
            schema_dir=SCHEMA_DIR,
        )

    except (ValueError, OSError, WorkflowError) as error:
        print(f"KASPberry input error: {error}", file=sys.stderr)
        return 2

    return run_snakemake(
        mode=args.mode,
        config_path=config_path,
        extra_args=snakemake_args,
    )


if __name__ == "__main__":
    sys.exit(main())
