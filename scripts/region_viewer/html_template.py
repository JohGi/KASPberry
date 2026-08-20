#!/usr/bin/env python3
# Author: Johanna Girodolle

"""HTML template rendering for the region overview viewer."""

from __future__ import annotations

import json
from pathlib import Path

from .constants import RESIZER_WIDTH, SIDEBAR_MIN_WIDTH, SIDEBAR_WIDTH
from .payload import build_config_payload

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
    config = build_config_payload()

    template = (_TEMPLATES_DIR / "region_viewer.html").read_text(encoding="utf-8")
    css = (_STATIC_DIR / "region_viewer.css").read_text(encoding="utf-8")
    js = "".join(
        (_STATIC_DIR / filename).read_text(encoding="utf-8")
        for filename in _REGION_VIEWER_JS_FILES
    )

    css = css.replace("{{ VIEWER_TOP_UI_HEIGHT }}", str(config["viewerTopUiHeight"]))
    css = css.replace("{{ SIDEBAR_WIDTH }}", str(SIDEBAR_WIDTH))
    css = css.replace("{{ SIDEBAR_MIN_WIDTH }}", str(SIDEBAR_MIN_WIDTH))
    css = css.replace("{{ RESIZER_WIDTH }}", str(RESIZER_WIDTH))

    js = js.replace("{{ REGION_DATA }}", json.dumps(region_data))
    js = js.replace("{{ CONFIG }}", json.dumps(config))

    html = template.replace("{{ CSS }}", css)
    html = html.replace("{{ JS }}", js)
    html = html.replace("{{ TITLE }}", region_viewer_title)
    html = html.replace("{{ VIEWER_TOP_UI_HEIGHT }}", str(config["viewerTopUiHeight"]))

    return html
