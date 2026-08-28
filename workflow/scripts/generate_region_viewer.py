#!/usr/bin/env python3
# Author: Johanna Girodolle

"""Generate an interactive Konva-based region overview HTML."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

from region_viewer.builder import RegionOverviewBuilder


LOGGER = logging.getLogger(__name__)


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate an interactive Konva-based region overview HTML."
    )
    parser.add_argument("--samples-tsv", type=Path, required=True)
    parser.add_argument("--mode", required=True, choices=["snps", "kasp"])
    parser.add_argument("--snp-summary", type=Path, required=True)
    parser.add_argument("--assay-summary", type=Path)
    parser.add_argument("--block-coords-tsv", type=Path, required=True)
    parser.add_argument("--snp-long", type=Path, required=True)
    parser.add_argument("--fasta-dir", type=Path, required=True)
    parser.add_argument("--summary-stats-json", type=Path, required=True)
    parser.add_argument("--mash-matrix", type=Path, required=True)
    parser.add_argument("--kimura2p-distmat-dir", type=Path, required=True)
    parser.add_argument("--masked-align-dir", type=Path, required=True)
    parser.add_argument("--masked-block-n-stats", type=Path, required=True)
    parser.add_argument("--gff-tracks-json", type=Path, required=True)
    parser.add_argument("--dotplot-manifest-json", type=Path, required=False)
    parser.add_argument("--analysis-settings-json", type=Path, required=True)
    parser.add_argument("--title", type=str, default="Region viewer")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
    )
    return parser.parse_args()


def validate_args(args: argparse.Namespace) -> None:
    """Validate mode-specific viewer inputs."""
    if args.mode == "kasp" and args.assay_summary is None:
        raise ValueError("--assay-summary is required in kasp mode.")
    

def setup_logging(log_level: str) -> None:
    """Configure logging."""
    logging.basicConfig(
        level=getattr(logging, log_level),
        format="%(levelname)s | %(message)s",
    )


def main() -> None:
    """Run the region overview generator."""
    args = parse_args()
    validate_args(args)
    setup_logging(args.log_level)

    builder = RegionOverviewBuilder(
        mode=args.mode,
        samples_tsv_path=args.samples_tsv,
        block_coords_tsv_path=args.block_coords_tsv,
        snp_long_path=args.snp_long,
        snp_summary_path=args.snp_summary,
        assay_summary_path=args.assay_summary,
        fasta_dir=args.fasta_dir,
        summary_stats_json_path=args.summary_stats_json,
        mash_matrix_path=args.mash_matrix,
        kimura2p_distmat_dir=args.kimura2p_distmat_dir,
        masked_align_dir=args.masked_align_dir,
        masked_block_n_stats_path=args.masked_block_n_stats,
        gff_tracks_json_path=args.gff_tracks_json,
        dotplot_manifest_json_path=args.dotplot_manifest_json,
        analysis_settings_json_path=args.analysis_settings_json,
        title=args.title,
        output_path=args.output,
    )
    builder.run()


if __name__ == "__main__":
    main()
