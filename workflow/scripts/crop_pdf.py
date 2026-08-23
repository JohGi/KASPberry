#!/usr/bin/env python3
# Author: Johanna Girodolle

"""Crop PDF pages using Ghostscript bounding boxes and PyMuPDF crop boxes."""

from __future__ import annotations

import argparse
import logging
import re
import subprocess
from pathlib import Path

import fitz


LOGGER = logging.getLogger(__name__)

HIRES_BBOX_PATTERN = re.compile(
    r"^%%HiResBoundingBox:\s+"
    r"(?P<x0>-?\d+(?:\.\d+)?)\s+"
    r"(?P<y0>-?\d+(?:\.\d+)?)\s+"
    r"(?P<x1>-?\d+(?:\.\d+)?)\s+"
    r"(?P<y1>-?\d+(?:\.\d+)?)"
)


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Crop PDF pages using Ghostscript bounding boxes."
    )
    parser.add_argument("--input", type=Path, required=True, help="Input PDF file.")
    parser.add_argument("--output", type=Path, required=True, help="Output cropped PDF file.")
    parser.add_argument(
        "--margin",
        type=float,
        default=0.0,
        help="Optional crop margin in PDF points. Default: 0.",
    )
    return parser.parse_args()


def configure_logging() -> None:
    """Configure application logging."""
    logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")


def run_ghostscript_bbox(input_path: Path) -> list[str]:
    """Run Ghostscript bbox device and return stderr lines."""
    command = [
        "gs",
        "-q",
        "-dNOPAUSE",
        "-dBATCH",
        "-sDEVICE=bbox",
        str(input_path),
    ]
    result = subprocess.run(
        command,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    return result.stderr.splitlines()


def parse_hires_bounding_boxes(lines: list[str]) -> list[fitz.Rect]:
    """Parse Ghostscript HiResBoundingBox lines."""
    boxes: list[fitz.Rect] = []

    for line in lines:
        match = HIRES_BBOX_PATTERN.match(line.strip())
        if match is None:
            continue

        boxes.append(
            fitz.Rect(
                float(match.group("x0")),
                float(match.group("y0")),
                float(match.group("x1")),
                float(match.group("y1")),
            )
        )

    return boxes


def convert_gs_bbox_to_pdf_rect(
    bbox: fitz.Rect,
    page_rect: fitz.Rect,
    margin: float,
) -> fitz.Rect:
    """Convert a Ghostscript bbox to a PyMuPDF crop rectangle."""
    rect = fitz.Rect(
        bbox.x0,
        page_rect.height - bbox.y1,
        bbox.x1,
        page_rect.height - bbox.y0,
    )

    return fitz.Rect(
        max(page_rect.x0, rect.x0 - margin),
        max(page_rect.y0, rect.y0 - margin),
        min(page_rect.x1, rect.x1 + margin),
        min(page_rect.y1, rect.y1 + margin),
    )


def crop_pdf(input_path: Path, output_path: Path, margin: float) -> None:
    """Crop all pages of a PDF using Ghostscript HiResBoundingBox values."""
    if not input_path.is_file():
        raise FileNotFoundError(f"Input PDF file not found: {input_path}")

    gs_lines = run_ghostscript_bbox(input_path)
    bounding_boxes = parse_hires_bounding_boxes(gs_lines)

    if not bounding_boxes:
        raise ValueError(f"No Ghostscript HiResBoundingBox found for {input_path}")

    output_path.parent.mkdir(parents=True, exist_ok=True)

    with fitz.open(input_path) as document:
        if len(bounding_boxes) != document.page_count:
            raise ValueError(
                "Ghostscript bounding box count does not match PDF page count: "
                f"{len(bounding_boxes)} bbox values for {document.page_count} pages."
            )

        for page, bbox in zip(document, bounding_boxes):
            crop_rect = convert_gs_bbox_to_pdf_rect(
                bbox=bbox,
                page_rect=page.rect,
                margin=margin,
            )
            page.set_cropbox(crop_rect)

        document.save(output_path)


def main() -> None:
    """Run PDF cropping."""
    configure_logging()
    args = parse_args()
    crop_pdf(
        input_path=args.input,
        output_path=args.output,
        margin=args.margin,
    )


if __name__ == "__main__":
    main()
