#!/usr/bin/env python3
"""Record provenance information for KASPberry runs."""

from __future__ import annotations

import copy
import importlib.metadata
import platform
import shlex
import socket
import subprocess
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import yaml


@dataclass
class RunRecord:
    """Files and metadata associated with one KASPberry run."""

    history_path: Path
    latest_path: Path
    data: dict
    started_at: datetime


def is_dry_run(snakemake_args: list[str]) -> bool:
    """Return True if Snakemake is being invoked in dry-run mode."""
    return any(
        arg in {"-n", "--dry-run", "--dryrun"}
        for arg in snakemake_args
    )


def _get_kaspberry_version(repo_root: Path) -> str:
    """
    Return the installed KASPberry version.

    Prefer package metadata when KASPberry is installed as a Python package.
    During development, fall back to the repository VERSION file.
    """
    try:
        return importlib.metadata.version("kaspberry")
    except importlib.metadata.PackageNotFoundError:
        pass

    version_file = repo_root / "VERSION"
    if version_file.is_file():
        version = version_file.read_text(encoding="utf-8").strip()
        if version:
            return version

    return "unknown"


def _get_snakemake_version() -> str:
    """Return the installed Snakemake version."""
    try:
        return importlib.metadata.version("snakemake")
    except importlib.metadata.PackageNotFoundError:
        return "unknown"


def _get_git_info(repo_root: Path) -> tuple[str | None, bool | None]:
    """
    Return the current Git commit and whether the working tree is dirty.

    Git information is optional and is unavailable for installations
    that are not Git working copies.
    """
    if not (repo_root / ".git").exists():
        return None, None

    try:
        commit = subprocess.run(
            [
                "git",
                "-C",
                str(repo_root),
                "rev-parse",
                "HEAD",
            ],
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()

        status = subprocess.run(
            [
                "git",
                "-C",
                str(repo_root),
                "status",
                "--porcelain",
            ],
            check=True,
            capture_output=True,
            text=True,
        ).stdout

    except (OSError, subprocess.CalledProcessError):
        return None, None

    return commit or None, bool(status.strip())


def _write_yaml(path: Path, data: dict) -> None:
    """Write YAML atomically."""
    path.parent.mkdir(parents=True, exist_ok=True)

    temporary_path = path.with_suffix(path.suffix + ".tmp")

    with temporary_path.open("w", encoding="utf-8") as handle:
        yaml.safe_dump(
            data,
            handle,
            sort_keys=False,
            allow_unicode=True,
        )

    temporary_path.replace(path)


def start_run_record(
    *,
    config: dict,
    mode: str,
    config_path: Path,
    repo_root: Path,
    argv: list[str],
    snakemake_args: list[str],
) -> RunRecord:
    """Create provenance files for a KASPberry run."""

    started_at = datetime.now().astimezone()

    run_id = started_at.strftime(
        "%Y%m%dT%H%M%S%f%z"
    )

    output_dir = Path(
        config["project"]["output_dir"]
    ).expanduser().resolve()

    history_path = (
        output_dir
        / "run_history"
        / f"{run_id}_{mode}.yaml"
    )

    latest_path = output_dir / "run_info.yaml"

    git_commit, git_dirty = _get_git_info(repo_root)

    data = {
        "run": {
            "id": run_id,
            "mode": mode,
            "status": "running",
            "started_at": started_at.isoformat(timespec="seconds"),
            "finished_at": None,
            "duration_seconds": None,
            "return_code": None,
        },
        "software": {
            "kaspberry_version": _get_kaspberry_version(repo_root),
            "git_commit": git_commit,
            "git_dirty": git_dirty,
            "snakemake_version": _get_snakemake_version(),
            "python_version": platform.python_version(),
        },
        "execution": {
            "command": shlex.join(["kaspberry", *argv]),
            "working_directory": str(Path.cwd()),
            "hostname": socket.gethostname(),
            "platform": platform.platform(),
            "snakemake_args": list(snakemake_args),
        },
        "inputs": {
            "config_file": str(config_path),
        },
        "resolved_config": copy.deepcopy(config),
    }

    record = RunRecord(
        history_path=history_path,
        latest_path=latest_path,
        data=data,
        started_at=started_at,
    )

    _write_yaml(record.history_path, record.data)
    _write_yaml(record.latest_path, record.data)

    return record


def finish_run_record(
    record: RunRecord,
    *,
    return_code: int,
    status: str | None = None,
) -> None:
    """Finalize provenance information after Snakemake exits."""

    finished_at = datetime.now().astimezone()

    if status is None:
        status = "success" if return_code == 0 else "failed"

    duration = (
        finished_at - record.started_at
    ).total_seconds()

    record.data["run"].update(
        {
            "status": status,
            "finished_at": finished_at.isoformat(timespec="seconds"),
            "duration_seconds": round(duration, 3),
            "return_code": return_code,
        }
    )

    _write_yaml(record.history_path, record.data)
    _write_yaml(record.latest_path, record.data)
