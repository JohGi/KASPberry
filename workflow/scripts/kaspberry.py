#!/usr/bin/env python3
"""Command-line interface for KASPberry."""

from __future__ import annotations

import argparse
import subprocess
import sys
import os
from pathlib import Path




import yaml
from snakemake_interface_common.exceptions import WorkflowError

from run_info import (
    finish_run_record,
    get_snakemake_working_directory,
    is_dry_run,
    start_run_record,
)
from validate_inputs import validate_inputs


REPO_ROOT = Path(__file__).resolve().parents[2]
SNAKEFILE = REPO_ROOT / "workflow" / "Snakefile"
SCHEMA_DIR = REPO_ROOT / "workflow" / "schemas"
DEFAULT_CONDA_PREFIX = REPO_ROOT / ".snakemake" / "conda"
DEFAULT_APPTAINER_PREFIX = REPO_ROOT / ".snakemake" / "apptainer"
DEFAULT_SOFTWARE_DEPLOYMENT_METHODS = ["conda", "apptainer"]


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


def get_conda_prefix_args(extra_args: list[str]) -> list[str]:
    """Return a shared Conda prefix unless the user selected one."""
    has_conda_prefix = any(
        arg == "--conda-prefix" or arg.startswith("--conda-prefix=")
        for arg in extra_args
    )
    if has_conda_prefix or "SNAKEMAKE_CONDA_PREFIX" in os.environ:
        return []

    return ["--conda-prefix", str(DEFAULT_CONDA_PREFIX)]


def get_apptainer_prefix_args(extra_args: list[str]) -> list[str]:
    """Return a shared Apptainer cache prefix unless the user selected one."""
    has_apptainer_prefix = any(
        arg in ("--apptainer-prefix", "--singularity-prefix")
        or arg.startswith(("--apptainer-prefix=", "--singularity-prefix="))
        for arg in extra_args
    )
    if has_apptainer_prefix or "APPTAINER_CACHEDIR" in os.environ:
        return []

    return ["--apptainer-prefix", str(DEFAULT_APPTAINER_PREFIX)]


def get_apptainer_args(extra_args: list[str]) -> list[str]:
    """Add a repository bind to forwarded Apptainer arguments.

    Snakemake runs from the analysis directory, while workflow shell commands
    can reference scripts using absolute paths below ``REPO_ROOT``.  The bind
    keeps those paths valid inside every Apptainer container.
    """
    repository_bind = f"{REPO_ROOT}:{REPO_ROOT}"
    bind_argument = f"--bind {repository_bind}"
    apptainer_args_indexes: list[tuple[int, bool]] = []

    for index, arg in enumerate(extra_args):
        if arg == "--apptainer-args":
            apptainer_args_indexes.append((index + 1, False))
        elif arg.startswith("--apptainer-args="):
            apptainer_args_indexes.append((index, True))

    supplied_apptainer_args = [
        (
            extra_args[index].split("=", 1)[1]
            if is_equals_form
            else extra_args[index]
        )
        for index, is_equals_form in apptainer_args_indexes
        if index < len(extra_args)
    ]
    if any(repository_bind in args for args in supplied_apptainer_args):
        return extra_args

    if not apptainer_args_indexes:
        return [*extra_args, "--apptainer-args", bind_argument]

    updated_args = extra_args.copy()
    for index, is_equals_form in apptainer_args_indexes:
        if index >= len(updated_args):
            continue

        existing_args = (
            updated_args[index].split("=", 1)[1]
            if is_equals_form
            else updated_args[index]
        )
        combined_args = " ".join(part for part in (existing_args, bind_argument) if part)
        updated_args[index] = (
            f"--apptainer-args={combined_args}"
            if is_equals_form
            else combined_args
        )

    return updated_args


def get_software_deployment_method_args(extra_args: list[str]) -> list[str]:
    """Enable KASPberry's Conda and Apptainer rule environments by default."""
    deployment_options = (
        "--software-deployment-method",
        "--deployment-method",
        "--deployment",
        "--sdm",
    )
    has_deployment_method = any(
        arg in deployment_options
        or arg.startswith(tuple(f"{option}=" for option in deployment_options))
        for arg in extra_args
    )
    if has_deployment_method:
        return []

    return ["--software-deployment-method", *DEFAULT_SOFTWARE_DEPLOYMENT_METHODS]


def run_snakemake(
    mode: str,
    config_path: Path,
    extra_args: list[str],
) -> int:
    """Run the requested Snakemake target."""
    apptainer_args = get_apptainer_args(extra_args)
    command = [
        "snakemake",
        "--snakefile",
        str(SNAKEFILE),
        "--configfile",
        str(config_path),
        *get_conda_prefix_args(extra_args),
        *get_apptainer_prefix_args(extra_args),
        *get_software_deployment_method_args(extra_args),
        *apptainer_args,
        mode,
    ]

    env = os.environ.copy()
    env["KASPBERRY_MODE"] = mode

    return subprocess.run(command, env=env).returncode


def run_kaspberry(
    *,
    mode: str,
    config: dict,
    config_path: Path,
    snakemake_args: list[str],
    working_directory: Path,
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
        working_directory=working_directory,
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

        working_directory = get_snakemake_working_directory(
            snakemake_args
        )

        validate_inputs(
            config=config,
            mode=args.mode,
            schema_dir=SCHEMA_DIR,
            working_directory=working_directory,
        )

        return run_kaspberry(
            mode=args.mode,
            config=config,
            config_path=config_path,
            snakemake_args=snakemake_args,
            working_directory=working_directory,
        )

    except (ValueError, OSError, WorkflowError) as error:
        print(f"KASPberry error: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
