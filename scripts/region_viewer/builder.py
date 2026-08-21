#!/usr/bin/env python3
# Author: Johanna Girodolle

"""Build the region overview HTML from workflow outputs."""

from __future__ import annotations

from pathlib import Path

from attrs import define, field

from .html_template import build_html
from .io import (
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
from .models import BlockFeature
from .payload import build_region_payload, build_sample_data
from .result_io import (
    group_assays_by_snp,
    read_assay_results,
    read_snp_results,
)
from .settings import read_analysis_settings


@define
class RegionOverviewBuilder:
    """Build the final Konva HTML from workflow outputs."""

    mode: str
    samples_tsv_path: Path
    block_coords_tsv_path: Path
    snp_long_path: Path
    snp_summary_path: Path
    fasta_dir: Path
    summary_stats_json_path: Path
    mash_matrix_path: Path
    kimura2p_distmat_dir: Path
    masked_block_n_stats_path: Path
    masked_align_dir: Path
    gff_tracks_json_path: Path
    title: str
    output_path: Path

    assay_summary_path: Path | None = field(default=None)
    config_yaml_path: Path | None = field(default=None)
    dotplot_manifest_json_path: Path | None = field(default=None)

    def run(self) -> None:
        """Run the full HTML generation workflow."""
        sample_records = read_samples(self.samples_tsv_path)
        sample_order = [record.sample for record in sample_records]

        fasta_lengths = read_fasta_lengths(self.fasta_dir)
        blocks_by_sample = read_blocks(self.block_coords_tsv_path)

        block_alignments = read_block_alignments(
            align_dir=self.masked_align_dir,
            block_ids=self.get_block_ids(blocks_by_sample),
        )

        snp_long = read_snp_long(self.snp_long_path)
        snps_by_sample = parse_snps(snp_long)

        snp_results = read_snp_results(
            path=self.snp_summary_path,
            mode=self.mode,
        )
        self.validate_snp_results(snps_by_sample, snp_results)

        assays_by_snp = {}

        if self.mode == "kasp":
            if self.assay_summary_path is None:
                raise ValueError(
                    "Assay summary path is required in KASP mode."
                )

            assays_by_snp = group_assays_by_snp(
                read_assay_results(self.assay_summary_path)
            )

        sample_data = build_sample_data(
            sample_records=sample_records,
            fasta_lengths=fasta_lengths,
            blocks_by_sample=blocks_by_sample,
            snps_by_sample=snps_by_sample,
        )

        gff_tracks_config = read_gff_tracks_json(
            self.gff_tracks_json_path
        )
        gff_tracks_by_sample = read_gff_gene_tracks(
            gff_tracks=gff_tracks_config,
            sample_data=sample_data,
        )

        region_data = build_region_payload(
            mode=self.mode,
            sample_data=sample_data,
            snp_results=snp_results,
            assays_by_snp=assays_by_snp,
            summary_stats=read_summary_stats(
                self.summary_stats_json_path
            ),
            mash_matrix=parse_mash_matrix(
                path=self.mash_matrix_path,
                sample_order=sample_order,
            ).to_dict(),
            kimura2p_matrices=parse_kimura2p_distmat_dir(
                distmat_dir=self.kimura2p_distmat_dir,
                sample_order=sample_order,
            ),
            masked_block_n_stats=read_masked_block_n_stats(
                self.masked_block_n_stats_path
            ),
            block_alignments=block_alignments,
            gff_tracks_by_sample=gff_tracks_by_sample,
            analysis_settings=read_analysis_settings(
                self.config_yaml_path,
                mode=self.mode,
            ),
            dotplots=read_dotplot_manifest(
                path=self.dotplot_manifest_json_path,
                output_path=self.output_path,
            ),
        )

        write_html(
            build_html(region_data, self.title),
            self.output_path,
        )

    @staticmethod
    def validate_snp_results(
        snps_by_sample,
        snp_results,
    ) -> None:
        """Ensure detected SNPs and aggregated results match."""
        displayed_ids = {
            snp.feature_id
            for sample_snps in snps_by_sample.values()
            for snp in sample_snps
        }
        result_ids = set(snp_results)

        missing = displayed_ids - result_ids
        unexpected = result_ids - displayed_ids

        if missing:
            raise ValueError(
                "SNPs missing from aggregated summary: "
                + ", ".join(sorted(missing)[:10])
            )

        if unexpected:
            raise ValueError(
                "Aggregated SNPs absent from snp_positions_long.tsv: "
                + ", ".join(sorted(unexpected)[:10])
            )

    @staticmethod
    def get_block_ids(
        blocks_by_sample: dict[str, list[BlockFeature]],
    ) -> list[str]:
        """Return sorted unique block IDs."""
        return sorted(
            {
                block.block_id
                for blocks in blocks_by_sample.values()
                for block in blocks
            }
        )
