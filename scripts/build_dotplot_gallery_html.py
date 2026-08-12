#!/usr/bin/env python3
# Author: Johanna Girodolle

"""Build a combined HTML gallery for pairwise dotplot SVG files."""

from __future__ import annotations

import argparse
import csv
import html
import logging
import os
from pathlib import Path

from attrs import define


LOGGER = logging.getLogger(__name__)


@define(frozen=True)
class SampleRecord:
    """Store one sample record from the sample sheet."""

    fasta: str
    sample: str


@define(frozen=True)
class GalleryCell:
    """Store one gallery cell in the matrix."""

    cell_type: str
    label: str = ""
    svg_rel_path: str = ""


@define(frozen=True)
class MatrixRow:
    """Store one matrix row with its vertical label."""

    row_label: str
    cells: list[GalleryCell]


@define(frozen=True)
class GalleryConfig:
    """Store gallery rendering configuration."""

    samples_path: Path
    svg_dir: Path
    output_path: Path
    pivot: str = ""
    title: str = ""

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Build a combined HTML gallery for pairwise dotplot SVG files."
    )
    parser.add_argument(
        "--samples",
        required=True,
        help="Input sample sheet with 2 or 3 tab-separated columns.",
    )
    parser.add_argument(
        "--svg-dir",
        required=True,
        help="Directory containing pairwise dotplot SVG files.",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output HTML file path.",
    )
    parser.add_argument(
        "--pivot",
        default="",
        help="Optional pivot sample name. If provided, render a pivot layout.",
    )
    parser.add_argument(
        "--title",
        default="",
        help="Optional title for the gallery.",
    )
    return parser.parse_args()


def configure_logging() -> None:
    """Configure application logging."""
    logging.basicConfig(
        level=logging.INFO,
        format="[%(levelname)s] %(message)s",
    )


def build_config_from_args(args: argparse.Namespace) -> GalleryConfig:
    """Build the gallery configuration."""
    return GalleryConfig(
        samples_path=Path(args.samples),
        svg_dir=Path(args.svg_dir),
        output_path=Path(args.output),
        pivot=str(args.pivot).strip(),
        title=str(args.title).strip(),
    )


def read_samples(samples_path: Path) -> list[SampleRecord]:
    """Read samples from the KASPberry genotype table."""
    records: list[SampleRecord] = []

    with open(samples_path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        if reader.fieldnames is None:
            raise ValueError(f"Genotype table has no header: {samples_path}")

        required_columns = {"genotype", "region_fasta"}
        missing_columns = required_columns - set(reader.fieldnames)

        if missing_columns:
            raise ValueError(
                f"Missing required columns in genotype table {samples_path}: "
                f"{sorted(missing_columns)}"
            )

        for row in reader:
            records.append(
                SampleRecord(
                    fasta=row["region_fasta"],
                    sample=row["genotype"],
                )
            )

    if not records:
        raise ValueError(f"Genotype table is empty: {samples_path}")

    return records


def extract_sample_names(records: list[SampleRecord]) -> list[str]:
    """Extract sample names while preserving input order."""
    return [record.sample for record in records]


def validate_pivot(sample_names: list[str], pivot: str) -> None:
    """Validate the optional pivot sample."""
    if pivot and pivot not in sample_names:
        raise ValueError(
            f"Unknown pivot sample {pivot!r}. Expected one of: {sample_names}"
        )


def build_pair_id(sample_a: str, sample_b: str) -> str:
    """Build a stable pair identifier."""
    return f"{sample_a}__vs__{sample_b}"


def build_svg_path(pair_id: str, svg_dir: Path) -> Path:
    """Build the expected SVG path for one pair."""
    return svg_dir / f"{pair_id}.svg"


def build_relative_svg_path(output_path: Path, svg_path: Path) -> str:
    """Build a relative SVG path from the output HTML file."""
    return os.path.relpath(svg_path, start=output_path.parent)


def create_empty_cell() -> GalleryCell:
    """Create an empty matrix cell."""
    return GalleryCell(cell_type="empty")


def create_plot_cell(svg_rel_path: str) -> GalleryCell:
    """Create a plot cell."""
    return GalleryCell(cell_type="plot", svg_rel_path=svg_rel_path)


def create_missing_cell(label: str) -> GalleryCell:
    """Create a missing plot cell."""
    return GalleryCell(cell_type="missing", label=label)


def build_full_matrix_column_headers(sample_names: list[str]) -> list[str]:
    """Build visible headers for the compact triangular layout."""
    return sample_names[1:]


def build_full_matrix_rows(
    sample_names: list[str],
    svg_dir: Path,
    output_path: Path,
) -> list[MatrixRow]:
    """Build compact upper-triangular rows without diagonal cells.

    For samples A, B, C, D, E:

        visible headers: B C D E

        row A: A vs B, A vs C, A vs D, A vs E
        row B: empty, B vs C, B vs D, B vs E
        row C: empty, empty, C vs D, C vs E
        row D: empty, empty, empty, D vs E
    """
    if len(sample_names) < 2:
        raise ValueError("At least two samples are required to build a dotplot matrix.")

    rows: list[MatrixRow] = []
    row_samples = sample_names[:-1]
    col_samples = sample_names[1:]
    sample_to_index = {sample: index for index, sample in enumerate(sample_names)}

    for row_sample in row_samples:
        row_index = sample_to_index[row_sample]
        row_cells: list[GalleryCell] = []

        for col_sample in col_samples:
            col_index = sample_to_index[col_sample]

            if col_index <= row_index:
                row_cells.append(create_empty_cell())
                continue

            pair_id = build_pair_id(row_sample, col_sample)
            svg_path = build_svg_path(pair_id, svg_dir)

            if svg_path.exists():
                row_cells.append(
                    create_plot_cell(
                        svg_rel_path=build_relative_svg_path(output_path, svg_path)
                    )
                )
            else:
                row_cells.append(
                    create_missing_cell(label=f"Missing SVG for {pair_id}")
                )

        rows.append(MatrixRow(row_label=row_sample, cells=row_cells))

    return rows


def build_pivot_column_headers(
    sample_names: list[str],
    pivot: str,
) -> list[str]:
    """Build visible headers for the pivot layout."""
    return [sample_name for sample_name in sample_names if sample_name != pivot]


def build_pivot_rows(
    sample_names: list[str],
    pivot: str,
    svg_dir: Path,
    output_path: Path,
) -> list[MatrixRow]:
    """Build pivot rows with one vertical row label and one row of plots."""
    row_cells: list[GalleryCell] = []

    for sample_name in sample_names:
        if sample_name == pivot:
            continue

        pair_id = build_pair_id(pivot, sample_name)
        svg_path = build_svg_path(pair_id, svg_dir)

        if svg_path.exists():
            row_cells.append(
                create_plot_cell(
                    svg_rel_path=build_relative_svg_path(output_path, svg_path)
                )
            )
        else:
            row_cells.append(
                create_missing_cell(label=f"Missing SVG for {pair_id}")
            )

    return [MatrixRow(row_label=pivot, cells=row_cells)]


def render_column_headers(column_headers: list[str]) -> str:
    """Render matrix column headers."""
    parts = ['<div class="corner-cell"></div>']

    for header in column_headers:
        parts.append(f'<div class="col-header">{html.escape(header)}</div>')

    return "\n".join(parts)


def render_row_label(row_label: str) -> str:
    """Render one vertical row label."""
    return (
        '<div class="row-label-cell">'
        f'<span>{html.escape(row_label)}</span>'
        "</div>"
    )


def render_gallery_cell(cell: GalleryCell) -> str:
    """Render one matrix cell."""
    if cell.cell_type == "empty":
        return '<div class="gallery-cell empty-cell"></div>'

    if cell.cell_type == "missing":
        return (
            '<div class="gallery-cell missing-cell">'
            f'<div class="missing-note">{html.escape(cell.label)}</div>'
            "</div>"
        )

    return (
        '<div class="gallery-cell plot-cell">'
        '<div class="plot-frame">'
        f'<img src="{html.escape(cell.svg_rel_path)}" alt="" loading="lazy" />'
        "</div>"
        "</div>"
    )


def render_rows(rows: list[MatrixRow]) -> str:
    """Render all matrix rows."""
    rendered_cells: list[str] = []

    for row in rows:
        rendered_cells.append(render_row_label(row.row_label))
        for cell in row.cells:
            rendered_cells.append(render_gallery_cell(cell))

    return "\n".join(rendered_cells)


def infer_full_matrix_grid_column_count(sample_names: list[str]) -> int:
    """Infer the grid column count for the compact triangular layout."""
    return len(sample_names)


def infer_pivot_grid_column_count(sample_names: list[str]) -> int:
    """Infer the grid column count for the pivot layout."""
    return len(sample_names)


def build_summary_label(sample_names: list[str], pivot: str) -> str:
    """Build the toolbar summary label."""
    if pivot:
        compared_samples = len(sample_names) - 1
        return (
            f"Mode: pivot | Pivot: {pivot} | "
            f"Compared samples: {compared_samples}"
        )
    return f"Mode: full matrix | Samples: {len(sample_names)}"

def build_html_document(
    sample_names: list[str],
    pivot: str,
    column_headers: list[str],
    rows: list[MatrixRow],
    title: str = "",
) -> str:
    """Build the full HTML document."""
    if pivot:
        grid_columns = infer_pivot_grid_column_count(sample_names)
    else:
        grid_columns = infer_full_matrix_grid_column_count(sample_names)

    header_html = render_column_headers(column_headers)
    summary_label = build_summary_label(sample_names, pivot)
    title_html = f"<h1>{html.escape(title)}</h1>" if title else ""

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1.0" />
  <title>Combined dotplot gallery</title>
  <style>
    :root {{
      --bg: #f3f4f6;
      --panel: #ffffff;
      --header-bg: #f3f4f6;
      --border: #d1d5db;
      --text: #111827;
      --muted: #6b7280;
      --shadow: 0 6px 18px rgba(15, 23, 42, 0.06);
      --cell-radius: 16px;
      --plot-min-height: 320px;
      --plot-col-width: 260px;
      --row-label-width: 72px;
      --grid-gap: 14px;
      --col-header-font-size: 28px;
      --row-label-font-size: 30px;
      --missing-font-size: 18px;
    }}

    * {{
      box-sizing: border-box;
    }}

    body {{
      margin: 0;
      font-family: Arial, sans-serif;
      background: var(--bg);
      color: var(--text);
    }}

    .page {{
      padding: 12px;
    }}

    .toolbar {{
      position: sticky;
      top: 0;
      z-index: 20;
      display: flex;
      gap: 10px;
      align-items: center;
      flex-wrap: wrap;
      padding: 4px 0 14px 0;
      background: var(--bg);
    }}

    .toolbar button {{
      border: 1px solid var(--border);
      background: var(--panel);
      color: var(--text);
      border-radius: 10px;
      padding: 8px 12px;
      cursor: pointer;
      box-shadow: var(--shadow);
      font-size: 14px;
    }}

    .summary {{
      color: var(--muted);
      font-size: 14px;
      margin-left: 4px;
    }}

    .viewport {{
      overflow-x: auto;
      overflow-y: hidden;
      border: 1px solid var(--border);
      border-radius: 18px;
      background: #f7f7f8;
      box-shadow: var(--shadow);
      padding: 16px;
    }}

    .zoom-layer {{
      position: relative;
      transform-origin: top left;
      height: auto;
    }}

    .gallery-grid {{
      display: grid;
      grid-template-columns: var(--row-label-width) repeat({grid_columns - 1}, minmax(var(--plot-col-width), 1fr));
      gap: var(--grid-gap);
      align-items: stretch;
      min-width: max-content;
      transform-origin: top left;
    }}

    .corner-cell,
    .col-header,
    .row-label-cell,
    .gallery-cell {{
      border-radius: var(--cell-radius);
    }}

    .corner-cell {{
      min-height: 50px;
      background: transparent;
      border: none;
      box-shadow: none;
    }}

    .col-header {{
      min-height: 50px;
      display: flex;
      align-items: center;
      justify-content: center;
      padding: 12px 14px;
      border: 1px solid var(--border);
      background: var(--header-bg);
      font-size: var(--col-header-font-size);
      font-weight: 700;
      text-align: center;
      line-height: 1;
    }}

    .row-label-cell {{
      min-height: var(--plot-min-height);
      height: 100%;
      display: flex;
      align-items: center;
      justify-content: center;
      padding: 10px 4px;
      border: 1px solid var(--border);
      background: var(--header-bg);
      box-shadow: 0 1px 0 rgba(17, 24, 39, 0.02);
      align-self: stretch;
    }}

    .row-label-cell span {{
      writing-mode: vertical-rl;
      transform: rotate(180deg);
      font-size: var(--row-label-font-size);
      font-weight: 700;
      line-height: 1;
      text-align: center;
      white-space: nowrap;
    }}

    .gallery-cell {{
      min-height: var(--plot-min-height);
      height: 100%;
    }}

    .plot-cell,
    .missing-cell {{
      border: 1px solid var(--border);
      background: var(--panel);
      box-shadow: 0 1px 0 rgba(17, 24, 39, 0.02);
    }}

    .empty-cell {{
      border: none;
      background: transparent;
      box-shadow: none;
    }}

    .missing-cell {{
      display: flex;
      align-items: center;
      justify-content: center;
      min-height: var(--plot-min-height);
      height: 100%;
      padding: 20px;
      text-align: center;
    }}

    .missing-note {{
      color: var(--muted);
      font-size: var(--missing-font-size);
      line-height: 1.25;
    }}

    .plot-cell {{
      padding: 10px;
    }}

    .plot-frame {{
      width: 100%;
      min-height: calc(var(--plot-min-height) - 20px);
      border-radius: 12px;
      background: #ffffff;
      display: flex;
      align-items: center;
      justify-content: center;
    }}

    .plot-frame img {{
      display: block;
      width: 100%;
      height: auto;
      user-select: none;
      pointer-events: none;
    }}

    h1 {{
      margin: 0 0 12px 0;
      font-size: 28px;
      font-weight: 700;
      color: var(--text);
    }}
  </style>
</head>
<body>
  <div class="page">
    {title_html}
    <div class="toolbar">
      <button type="button" id="zoom-out">-</button>
      <button type="button" id="zoom-in">+</button>
      <button type="button" id="zoom-reset">Reset</button>
      <span class="summary">{html.escape(summary_label)}</span>
    </div>

    <div class="viewport">
      <div class="zoom-layer" id="zoom-layer">
        <div class="gallery-grid" id="gallery-grid">
          {header_html}
          {render_rows(rows)}
        </div>
      </div>
    </div>
  </div>

  <script>
    const zoomLayer = document.getElementById("zoom-layer");
    const galleryGrid = document.getElementById("gallery-grid");
    const zoomInButton = document.getElementById("zoom-in");
    const zoomOutButton = document.getElementById("zoom-out");
    const zoomResetButton = document.getElementById("zoom-reset");

    const defaultZoomLevel = 0.4;
    let zoomLevel = defaultZoomLevel;
    const zoomStep = 0.1;
    const zoomMin = 0.1;
    const zoomMax = 4.0;

    function applyZoom() {{
      galleryGrid.style.transform = `scale(${{zoomLevel}})`;

      const rawWidth = galleryGrid.scrollWidth;
      const rawHeight = galleryGrid.scrollHeight;

      zoomLayer.style.width = `${{rawWidth * zoomLevel}}px`;
      zoomLayer.style.height = `${{rawHeight * zoomLevel}}px`;
    }}

    function zoomIn() {{
      zoomLevel = Math.min(zoomMax, zoomLevel + zoomStep);
      applyZoom();
    }}

    function zoomOut() {{
      zoomLevel = Math.max(zoomMin, zoomLevel - zoomStep);
      applyZoom();
    }}

    function zoomReset() {{
      zoomLevel = defaultZoomLevel;
      applyZoom();
    }}

    function applyZoomAfterLayout() {{
      requestAnimationFrame(applyZoom);
    }}

    zoomInButton.addEventListener("click", zoomIn);
    zoomOutButton.addEventListener("click", zoomOut);
    zoomResetButton.addEventListener("click", zoomReset);
    window.addEventListener("resize", applyZoomAfterLayout);
    window.addEventListener("load", applyZoomAfterLayout);

    document.querySelectorAll(".plot-frame img").forEach((image) => {{
      image.addEventListener("load", applyZoomAfterLayout);
    }});

    applyZoomAfterLayout();
  </script>
</body>
</html>
"""


def write_html(output_path: Path, content: str) -> None:
    """Write the HTML document to disk."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(content, encoding="utf-8")


def main() -> None:
    """Run the gallery builder."""
    configure_logging()
    args = parse_args()
    config = build_config_from_args(args)
    sample_records = read_samples(config.samples_path)
    sample_names = extract_sample_names(sample_records)
    validate_pivot(sample_names, config.pivot)

    if config.pivot:
        column_headers = build_pivot_column_headers(sample_names, config.pivot)
        rows = build_pivot_rows(
            sample_names=sample_names,
            pivot=config.pivot,
            svg_dir=config.svg_dir,
            output_path=config.output_path,
        )
    else:
        column_headers = build_full_matrix_column_headers(sample_names)
        rows = build_full_matrix_rows(
            sample_names=sample_names,
            svg_dir=config.svg_dir,
            output_path=config.output_path,
        )

    html_document = build_html_document(
        sample_names=sample_names,
        pivot=config.pivot,
        column_headers=column_headers,
        rows=rows,
        title=config.title,
    )
    write_html(config.output_path, html_document)
    LOGGER.info("Wrote gallery HTML to %s", config.output_path)


if __name__ == "__main__":
    main()
