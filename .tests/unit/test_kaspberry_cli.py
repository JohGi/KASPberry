from __future__ import annotations

import sys
from pathlib import Path
from types import SimpleNamespace

import pytest


REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "workflow" / "scripts"))

import kaspberry


def run_snakemake_and_capture_command(
    monkeypatch: pytest.MonkeyPatch,
    extra_args: list[str],
) -> list[str]:
    """Run the wrapper with a mocked Snakemake process."""
    captured_command: list[str] = []

    def mock_run(command: list[str], **_kwargs: object) -> SimpleNamespace:
        captured_command.extend(command)
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(kaspberry.subprocess, "run", mock_run)

    assert kaspberry.run_snakemake(
        mode="snps",
        config_path=Path("config.yaml"),
        extra_args=extra_args,
    ) == 0

    return captured_command


def test_run_snakemake_adds_shared_conda_prefix_by_default(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The wrapper shares Conda environments from the repository root."""
    monkeypatch.delenv("SNAKEMAKE_CONDA_PREFIX", raising=False)

    command = run_snakemake_and_capture_command(monkeypatch, ["--cores", "1"])

    assert command.count("--conda-prefix") == 1
    assert command[command.index("--conda-prefix") + 1] == str(
        REPO_ROOT / ".snakemake" / "conda"
    )


@pytest.mark.parametrize(
    "extra_args",
    [
        ["--conda-prefix", "/custom/conda"],
        ["--conda-prefix=/custom/conda"],
    ],
)
def test_run_snakemake_respects_explicit_conda_prefix(
    monkeypatch: pytest.MonkeyPatch,
    extra_args: list[str],
) -> None:
    """Forwarded Conda prefix options take precedence over the default."""
    monkeypatch.delenv("SNAKEMAKE_CONDA_PREFIX", raising=False)

    command = run_snakemake_and_capture_command(monkeypatch, extra_args)

    assert all(arg in command for arg in extra_args)
    assert command.count("--conda-prefix") == (
        1 if extra_args[0] == "--conda-prefix" else 0
    )
    assert str(REPO_ROOT / ".snakemake" / "conda") not in command


def test_run_snakemake_respects_conda_prefix_environment_variable(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """SNAKEMAKE_CONDA_PREFIX suppresses the wrapper default."""
    monkeypatch.setenv("SNAKEMAKE_CONDA_PREFIX", "/custom/conda")

    command = run_snakemake_and_capture_command(monkeypatch, ["--cores", "1"])

    assert "--conda-prefix" not in command
    assert str(REPO_ROOT / ".snakemake" / "conda") not in command


def test_run_snakemake_adds_shared_apptainer_prefix_by_default(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The wrapper shares Apptainer images from the repository root."""
    monkeypatch.delenv("APPTAINER_CACHEDIR", raising=False)

    command = run_snakemake_and_capture_command(monkeypatch, ["--cores", "1"])

    assert command.count("--apptainer-prefix") == 1
    assert command[command.index("--apptainer-prefix") + 1] == str(
        REPO_ROOT / ".snakemake" / "apptainer"
    )


@pytest.mark.parametrize(
    "extra_args",
    [
        ["--apptainer-prefix", "/custom/apptainer"],
        ["--apptainer-prefix=/custom/apptainer"],
        ["--singularity-prefix", "/custom/singularity"],
        ["--singularity-prefix=/custom/singularity"],
    ],
)
def test_run_snakemake_respects_explicit_apptainer_prefix(
    monkeypatch: pytest.MonkeyPatch,
    extra_args: list[str],
) -> None:
    """Forwarded Apptainer and Singularity prefixes override the default."""
    monkeypatch.delenv("APPTAINER_CACHEDIR", raising=False)

    command = run_snakemake_and_capture_command(monkeypatch, extra_args)

    assert all(arg in command for arg in extra_args)
    assert str(REPO_ROOT / ".snakemake" / "apptainer") not in command


def test_run_snakemake_respects_apptainer_cache_environment_variable(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """APPTAINER_CACHEDIR suppresses the wrapper default."""
    monkeypatch.setenv("APPTAINER_CACHEDIR", "/custom/cache")

    command = run_snakemake_and_capture_command(monkeypatch, ["--cores", "1"])

    assert "--apptainer-prefix" not in command
    assert str(REPO_ROOT / ".snakemake" / "apptainer") not in command


def test_run_snakemake_adds_default_software_deployment_methods(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Conda and Apptainer environments are enabled without extra user flags."""
    command = run_snakemake_and_capture_command(monkeypatch, ["--cores", "1"])

    index = command.index("--software-deployment-method")
    assert command[index + 1 : index + 3] == ["conda", "apptainer"]


@pytest.mark.parametrize(
    "extra_args",
    [
        ["--sdm", "conda"],
        ["--sdm=conda"],
        ["--software-deployment-method", "apptainer"],
        ["--software-deployment-method=apptainer"],
        ["--deployment-method", "conda"],
        ["--deployment=apptainer"],
    ],
)
def test_run_snakemake_respects_explicit_software_deployment_methods(
    monkeypatch: pytest.MonkeyPatch,
    extra_args: list[str],
) -> None:
    """A user-selected Snakemake deployment method is forwarded unchanged."""
    command = run_snakemake_and_capture_command(monkeypatch, extra_args)

    assert all(arg in command for arg in extra_args)
    assert command.count("--software-deployment-method") == (
        1 if extra_args[0] == "--software-deployment-method" else 0
    )


def test_run_snakemake_adds_repository_bind_by_default(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The repository is visible at the same path in containers."""
    command = run_snakemake_and_capture_command(monkeypatch, ["--cores", "1"])

    assert command[command.index("--apptainer-args") + 1] == (
        f"--bind {REPO_ROOT}:{REPO_ROOT}"
    )


@pytest.mark.parametrize(
    ("extra_args", "expected_args"),
    [
        (
            ["--apptainer-args", "--containall"],
            "--containall",
        ),
        (
            ["--apptainer-args=--containall"],
            "--containall",
        ),
    ],
)
def test_run_snakemake_appends_repository_bind_to_apptainer_args(
    monkeypatch: pytest.MonkeyPatch,
    extra_args: list[str],
    expected_args: str,
) -> None:
    """User Apptainer flags remain intact when the repository bind is added."""
    command = run_snakemake_and_capture_command(monkeypatch, extra_args)
    bind_argument = f"--bind {REPO_ROOT}:{REPO_ROOT}"

    if extra_args[0] == "--apptainer-args":
        apptainer_args = command[command.index("--apptainer-args") + 1]
    else:
        apptainer_args = next(
            arg.split("=", 1)[1]
            for arg in command
            if arg.startswith("--apptainer-args=")
        )

    assert apptainer_args == f"{expected_args} {bind_argument}"


def test_run_snakemake_does_not_duplicate_existing_repository_bind(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """An explicit repository bind is left unchanged."""
    existing_args = f"--containall --bind {REPO_ROOT}:{REPO_ROOT}"
    command = run_snakemake_and_capture_command(
        monkeypatch,
        ["--apptainer-args", existing_args],
    )

    assert command[command.index("--apptainer-args") + 1] == existing_args


def test_run_snakemake_preserves_unrelated_extra_arguments(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Unrelated Snakemake arguments retain their original values and order."""
    extra_args = ["--cores", "4", "--directory", "/analysis", "--dry-run"]
    command = run_snakemake_and_capture_command(monkeypatch, extra_args)

    positions = [command.index(arg) for arg in extra_args]
    assert positions == sorted(positions)
    assert [command[position] for position in positions] == extra_args
