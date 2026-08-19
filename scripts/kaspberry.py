#!/usr/bin/env python3
"""Command-line interface for KASPberry."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import yaml
from snakemake_interface_common.exceptions import WorkflowError

from run_info import (
    finish_run_record,
    is_dry_run,
    start_run_record,
)
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


def run_kaspberry(
    *,
    mode: str,
    config: dict,
    config_path: Path,
    snakemake_args: list[str],
) -> int:
    """Run KASPberry and record provenance for non-dry runs."""
    if is_dry_run(snakemake_args):
        return run_snakemake(
            mode=mode,
            config_path=config_path,
            extra_args=snakemake_args,
        )

    run_record = start_run_record(
        config=config,
        mode=mode,
        config_path=config_path,
        repo_root=REPO_ROOT,
        argv=sys.argv[1:],
        snakemake_args=snakemake_args,
    )

    try:
        return_code = run_snakemake(
            mode=mode,
            config_path=config_path,
            extra_args=snakemake_args,
        )

    except KeyboardInterrupt:
        finish_run_record(
            run_record,
            return_code=130,
            status="interrupted",
        )
        return 130

    except OSError:
        finish_run_record(
            run_record,
            return_code=1,
            status="failed",
        )
        raise

    finish_run_record(
        run_record,
        return_code=return_code,
    )

    return return_code


def main() -> int:
    """Run the KASPberry CLI."""
    args, snakemake_args = parse_args()

    try:
        config_path = args.config.resolve()

        config = load_config(config_path)

        validate_inputs(
            config=config,
            mode=args.mode,
            schema_dir=SCHEMA_DIR,
        )

        return run_kaspberry(
            mode=args.mode,
            config=config,
            config_path=config_path,
            snakemake_args=snakemake_args,
        )

    except (ValueError, OSError, WorkflowError) as error:
        print(f"KASPberry error: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
