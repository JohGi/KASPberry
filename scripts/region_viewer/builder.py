#!/usr/bin/env python3
# Author: Johanna Girodolle

"""Build the region overview HTML from workflow outputs."""

from __future__ import annotations

import logging
from pathlib import Path

import yaml
from attrs import define, field

from .html_template import build_html
from .io import (
    count_unique_snps,
    parse_kimura2p_distmat_dir,
    parse_mash_matrix,
    parse_snps,
    read_block_alignments,
    read_blocks,
    read_dotplot_manifest,
    read_fasta_lengths,
    read_gff_gene_tracks,
    read_gff_tracks_json,
    read_masked_block_n_stats,
    read_samples,
    read_snp_long,
    read_summary_stats,
    write_html,
)

LOGGER = logging.getLogger(__name__)
from .models import BlockFeature
from .payload import build_region_payload, build_sample_data


import csv


def read_analysis_settings(config_yaml_path: Path | None) -> dict[str, object]:
    """Extract display-relevant analysis settings from the pipeline config."""
    if config_yaml_path is None:
        return {}

    try:
        raw = yaml.safe_load(config_yaml_path.read_text(encoding="utf-8")) or {}
    except Exception:
        LOGGER.warning("Could not read config YAML: %s", config_yaml_path)
        return {}

    snps = raw.get("snps") or {}

    min_len = snps.get("min_block_length")
    min_flank = snps.get("min_snp_flank")

    genotype_path = Path(raw["inputs"]["genotypes"])

    groups: dict[str, list[str]] = {}

    with genotype_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        for row in reader:
            genotype = row["genotype"].strip()
            group = row["group"].strip()

            if group:
                groups.setdefault(group, []).append(genotype)

    return {
        "minimum_block_length_bp": int(min_len) if min_len is not None else None,
        "minimum_snp_flank_bp": int(min_flank) if min_flank is not None else None,
        "snp_groups": groups,
    }

@define
class RegionOverviewBuilder:
    """Build the final Konva HTML from workflow outputs."""

    samples_tsv_path: Path
    block_coords_tsv_path: Path
    snp_long_path: Path
    fasta_dir: Path
    summary_stats_json_path: Path
    mash_matrix_path: Path
    kimura2p_distmat_dir: Path
    masked_block_n_stats_path: Path
    masked_align_dir: Path
    gff_tracks_json_path: Path
    title: str
    output_path: Path
    config_yaml_path: Path | None = field(default=None)
    dotplot_manifest_json_path: Path | None = field(default=None)

    def run(self) -> None:
        """Run the full HTML generation workflow."""
        sample_records = read_samples(self.samples_tsv_path)
        sample_order = [record.sample for record in sample_records]

        fasta_lengths = read_fasta_lengths(self.fasta_dir)
        blocks_by_sample = read_blocks(self.block_coords_tsv_path)
        block_ids = self.get_block_ids(blocks_by_sample)
        block_alignments = read_block_alignments(
            align_dir=self.masked_align_dir,
            block_ids=block_ids,
        )
        snp_long = read_snp_long(self.snp_long_path)
        snps_by_sample = parse_snps(snp_long)

        summary_stats = read_summary_stats(self.summary_stats_json_path)
        summary_stats.setdefault("global", {})["n_snps_kept"] = count_unique_snps(snp_long)
        mash_matrix = parse_mash_matrix(
            path=self.mash_matrix_path,
            sample_order=sample_order,
        ).to_dict()
        kimura2p_matrices = parse_kimura2p_distmat_dir(
            distmat_dir=self.kimura2p_distmat_dir,
            sample_order=sample_order,
        )
        masked_block_n_stats = read_masked_block_n_stats(
            self.masked_block_n_stats_path
        )

        sample_data = build_sample_data(
            sample_records=sample_records,
            fasta_lengths=fasta_lengths,
            blocks_by_sample=blocks_by_sample,
            snps_by_sample=snps_by_sample,
        )

        gff_tracks_config = read_gff_tracks_json(self.gff_tracks_json_path)
        gff_tracks_by_sample = read_gff_gene_tracks(
            gff_tracks=gff_tracks_config,
            sample_data=sample_data,
        )

        analysis_settings = read_analysis_settings(self.config_yaml_path)

        dotplot_records = read_dotplot_manifest(
            path=self.dotplot_manifest_json_path,
            output_path=self.output_path,
        )

        region_data = build_region_payload(
            sample_data=sample_data,
            summary_stats=summary_stats,
            mash_matrix=mash_matrix,
            kimura2p_matrices=kimura2p_matrices,
            masked_block_n_stats=masked_block_n_stats,
            block_alignments=block_alignments,
            gff_tracks_by_sample=gff_tracks_by_sample,
            analysis_settings=analysis_settings,
            dotplots=dotplot_records,
        )

        html = build_html(region_data, self.title)
        write_html(html, self.output_path)

    @staticmethod
    def get_block_ids(blocks_by_sample: dict[str, list[BlockFeature]]) -> list[str]:
        """Return sorted unique block IDs from sample-indexed block records."""
        return sorted(
            {
                block.block_id
                for sample_blocks in blocks_by_sample.values()
                for block in sample_blocks
            }
        )
