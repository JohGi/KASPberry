#!/usr/bin/env python3
# Author: Johanna Girodolle

"""Read display-relevant KASPberry analysis settings."""

from __future__ import annotations

import json
from pathlib import Path

def read_analysis_settings(
    analysis_settings_json_path: Path,
) -> dict[str, object]:
    """Read the display settings resolved by the Snakemake workflow."""
    try:
        settings = json.loads(
            analysis_settings_json_path.read_text(encoding="utf-8")
        )
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(
            "Could not read analysis settings JSON: "
            f"{analysis_settings_json_path}"
        ) from error

    if not isinstance(settings, dict):
        raise ValueError(
            "Analysis settings JSON must contain an object: "
            f"{analysis_settings_json_path}"
        )

    return settings
