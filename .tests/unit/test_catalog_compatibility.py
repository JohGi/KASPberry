from __future__ import annotations

from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[2]
RUN_SUMMARY = REPO_ROOT / "workflow" / "rules" / "run_summary.smk"


def load_get_run_summary_mode(
    monkeypatch: pytest.MonkeyPatch,
    mode: str | None,
):
    """Load the rule helper with a controlled KASPBERRY_MODE value."""
    if mode is None:
        monkeypatch.delenv("KASPBERRY_MODE", raising=False)
    else:
        monkeypatch.setenv("KASPBERRY_MODE", mode)

    helpers, _, _ = RUN_SUMMARY.read_text(encoding="utf-8").partition(
        "\nrule write_run_summary:"
    )
    namespace: dict[str, object] = {}
    exec(helpers, namespace)
    return namespace["get_run_summary_mode"]


def test_run_summary_mode_defaults_to_snps(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    assert load_get_run_summary_mode(monkeypatch, None)() == "snps"


@pytest.mark.parametrize("mode", ["snps", "kasp"])
def test_run_summary_mode_preserves_cli_mode(
    monkeypatch: pytest.MonkeyPatch,
    mode: str,
) -> None:
    assert load_get_run_summary_mode(monkeypatch, mode)() == mode


def test_run_summary_mode_rejects_invalid_explicit_value(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    get_run_summary_mode = load_get_run_summary_mode(monkeypatch, "invalid")

    with pytest.raises(ValueError, match="Invalid KASPBERRY_MODE"):
        get_run_summary_mode()
