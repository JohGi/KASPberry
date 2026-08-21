#!/usr/bin/env python3
# Author: Johanna Girodolle

"""HTML template rendering for the region overview viewer."""

from __future__ import annotations

import json
from pathlib import Path

_TEMPLATES_DIR = Path(__file__).parent / "templates"
_STATIC_DIR = Path(__file__).parent / "static"
_REGION_VIEWER_JS_FILES = (
    "core.js",
    "sidebar.js",
    "region.js",
    "alignment.js",
    "controls.js",
    "dotplot.js",
    "main.js",
)


def build_html(region_data: dict[str, object], region_viewer_title: str = "Region viewer") -> str:
    """Render the final HTML document from external template and static files."""
    template = (
        _TEMPLATES_DIR / "region_viewer.html"
    ).read_text(encoding="utf-8")
    css = (_STATIC_DIR / "region_viewer.css").read_text(encoding="utf-8")
    js = "\n".join(
        (_STATIC_DIR / filename).read_text(encoding="utf-8")
        for filename in _REGION_VIEWER_JS_FILES
    )

    js = js.replace("{{ REGION_DATA }}", json.dumps(region_data))

    html = template.replace("{{ CSS }}", css)
    html = html.replace("{{ JS }}", js)
    html = html.replace("{{ TITLE }}", region_viewer_title)
    return html
