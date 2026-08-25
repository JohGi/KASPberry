from __future__ import annotations

import filecmp
import shutil
import subprocess
from pathlib import Path

import yaml

def copy_case_input(case_input_dir: Path, workdir: Path) -> None:
    """Copy one integration test input case into a temporary work directory."""
    shutil.copytree(case_input_dir, workdir, dirs_exist_ok=True)

def copy_shared_resources(resources_dir: Path, workdir: Path) -> None:
    """Copy shared integration test resources into the temporary work directory."""
    if not resources_dir.exists():
        return

    for path in resources_dir.iterdir():
        destination = workdir / path.name
        if path.is_dir():
            shutil.copytree(path, destination, dirs_exist_ok=True)
        else:
            shutil.copy2(path, destination)

def deep_update(base: dict, override: dict) -> dict:
    """Recursively update a nested dictionary."""
    result = dict(base)
    for key, value in override.items():
        if (
            key in result
            and isinstance(result[key], dict)
            and isinstance(value, dict)
        ):
            result[key] = deep_update(result[key], value)
        else:
            result[key] = value
    return result


def prepare_test_config(
    base_config: Path,
    workdir: Path,
    override_config: Path | None = None,
) -> Path:
    """Create a test config file in the temporary work directory."""
    with base_config.open("r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle)

    if override_config is not None and override_config.exists():
        with override_config.open("r", encoding="utf-8") as handle:
            override = yaml.safe_load(handle)
        config = deep_update(config, override)

    config_dir = workdir / "config"
    config_dir.mkdir(parents=True, exist_ok=True)

    config_path = config_dir / "config.yaml"
    with config_path.open("w", encoding="utf-8") as handle:
        yaml.safe_dump(config, handle, sort_keys=False)

    return config_path

def run_snakemake(
    repo_root: Path,
    workdir: Path,
    target: str,
    configfile: Path,
) -> None:
    """Run Snakemake on one target in an isolated work directory."""
    subprocess.run(
        [
            "snakemake",
            "-s",
            str(repo_root / "workflow" / "Snakefile"),
            "--directory",
            str(workdir),
            "--configfile",
            str(configfile),
            "--cores",
            "1",
            "--sdm",
            "conda",
            "apptainer",
            "--conda-prefix",
            str(repo_root / ".snakemake" / "conda"),
            "--apptainer-prefix",
            str(repo_root / ".snakemake" / "apptainer"),
            "--apptainer-args",
            f"--bind {repo_root}:{repo_root}",
            target,
        ],
        check=True,
    )

def should_ignore(
    path: Path,
    ignore_names: set[str] | None = None,
    ignore_suffixes: set[str] | None = None,
) -> bool:
    """Return True if a file should be ignored during directory comparison."""
    if ignore_names and path.name in ignore_names:
        return True

    if ignore_suffixes and any(path.name.endswith(suffix) for suffix in ignore_suffixes):
        return True

    return False


def list_files(
    root_dir: Path,
    ignore_names: set[str] | None = None,
    ignore_suffixes: set[str] | None = None,
) -> list[Path]:
    """List all non-ignored files relative to a root directory."""
    return sorted(
        path.relative_to(root_dir)
        for path in root_dir.rglob("*")
        if path.is_file()
        and not should_ignore(
            path=path,
            ignore_names=ignore_names,
            ignore_suffixes=ignore_suffixes,
        )
    )


def compare_directories(
    expected_dir: Path,
    observed_dir: Path,
    ignore_names: set[str] | None = None,
    ignore_suffixes: set[str] | None = None,
) -> None:
    """Recursively compare two directories while optionally ignoring some files."""
    assert expected_dir.exists(), f"Missing expected directory: {expected_dir}"
    assert observed_dir.exists(), f"Missing observed directory: {observed_dir}"

    expected_files = list_files(
        root_dir=expected_dir,
        ignore_names=ignore_names,
        ignore_suffixes=ignore_suffixes,
    )
    observed_files = list_files(
        root_dir=observed_dir,
        ignore_names=ignore_names,
        ignore_suffixes=ignore_suffixes,
    )

    assert expected_files == observed_files, (
        f"File lists differ.\nExpected: {expected_files}\nObserved: {observed_files}"
    )

    for relative_path in expected_files:
        expected_file = expected_dir / relative_path
        observed_file = observed_dir / relative_path
        assert filecmp.cmp(expected_file, observed_file, shallow=False), (
            f"Files differ: {relative_path}"
        )
