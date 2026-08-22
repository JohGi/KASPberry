#!/usr/bin/env python3
# Author: Johanna Girodolle

from __future__ import annotations

from .models import (
    AssayResult,
    BlockAlignment,
    BlockFeature,
    DotplotRecord,
    GffTrack,
    SampleData,
    SampleRecord,
    SnpFeature,
    SnpResult,
)


def build_sample_data(
    sample_records: list[SampleRecord],
    fasta_lengths: dict[str, int],
    blocks_by_sample: dict[str, list[BlockFeature]],
    snps_by_sample: dict[str, list[SnpFeature]],
) -> list[SampleData]:
    """Build plotting data in samples TSV order."""
    sample_data = []

    for record in sample_records:
        sample = record.sample

        if sample not in fasta_lengths:
            raise ValueError(
                f"Missing FASTA length for sample {sample!r}"
            )

        sample_data.append(
            SampleData(
                sample=sample,
                region_length=fasta_lengths[sample],
                region_start_in_source_seq=record.region_start_in_source_seq,
                blocks=sorted(
                    blocks_by_sample.get(sample, []),
                    key=lambda block: (
                        block.block_start_in_region,
                        block.block_end_in_region,
                        block.block_id,
                    ),
                ),
                snps=sorted(
                    snps_by_sample.get(sample, []),
                    key=lambda snp: (
                        snp.pos_in_region,
                        snp.block_id,
                        snp.aln_pos,
                    ),
                ),
            )
        )

    return sample_data


def build_region_payload(
    sample_data: list[SampleData],
    mode: str = "snps",
    snp_results: dict[str, SnpResult] | None = None,
    assays_by_snp: dict[str, list[AssayResult]] | None = None,
    summary_stats: dict[str, object] | None = None,
    mash_matrix: dict[str, object] | None = None,
    kimura2p_matrices: dict[str, dict[str, object]] | None = None,
    masked_block_n_stats: dict | None = None,
    block_alignments: dict[str, BlockAlignment] | None = None,
    gff_tracks_by_sample: dict[str, list[GffTrack]] | None = None,
    analysis_settings: dict[str, object] | None = None,
    dotplots: list[DotplotRecord] | None = None,
) -> dict[str, object]:
    """Build the JSON payload injected into the HTML."""
    gff_tracks_by_sample = gff_tracks_by_sample or {}

    return {
        "title": "Region overview",
        "mode": mode,
        "max_region_length": max(
            sample.region_length
            for sample in sample_data
        ),
        "samples": [
            {
                "sample": sample.sample,
                "region_length": sample.region_length,
                "region_start_in_source_seq":
                    sample.region_start_in_source_seq,
                "blocks": [
                    {
                        "feature_id": block.feature_id,
                        "block_id": block.block_id,
                        "block_start_in_region":
                            block.block_start_in_region,
                        "block_end_in_region":
                            block.block_end_in_region,
                        "block_start_in_source_seq":
                            block.block_start_in_source_seq,
                        "block_end_in_source_seq":
                            block.block_end_in_source_seq,
                    }
                    for block in sample.blocks
                ],
                "snps": [
                    {
                        "feature_id": snp.feature_id,
                        "block_id": snp.block_id,
                        "aln_pos": snp.aln_pos,
                        "nt": snp.nt,
                        "pos_in_block": snp.pos_in_block,
                        "pos_in_region": snp.pos_in_region,
                        "pos_in_source_seq": snp.pos_in_source_seq,
                    }
                    for snp in sample.snps
                ],
                "gff_tracks": [
                    {
                        "sample": track.sample,
                        "track_name": track.track_name,
                        "features": [
                            {
                                "sample": feature.sample,
                                "track_name": feature.track_name,
                                "gene_id": feature.gene_id,
                                "source_seq_id": feature.source_seq_id,
                                "start_in_source_seq":
                                    feature.start_in_source_seq,
                                "end_in_source_seq":
                                    feature.end_in_source_seq,
                                "start_in_region":
                                    feature.start_in_region,
                                "end_in_region":
                                    feature.end_in_region,
                                "attributes": feature.attributes,
                                "strand": feature.strand,
                            }
                            for feature in track.features
                        ],
                    }
                    for track in gff_tracks_by_sample.get(
                        sample.sample, []
                    )
                ],
            }
            for sample in sample_data
        ],
        "snp_results": {
            snp_id: result.to_payload()
            for snp_id, result in (snp_results or {}).items()
        },
        "assays_by_snp": {
            snp_id: [
                assay.to_payload()
                for assay in assays
            ]
            for snp_id, assays in (assays_by_snp or {}).items()
        },
        "summary_stats": summary_stats or {},
        "mash_matrix": mash_matrix or {},
        "kimura2p_matrices": kimura2p_matrices or {},
        "masked_block_n_stats": masked_block_n_stats or {},
        "block_alignments": {
            block_id: alignment.to_payload()
            for block_id, alignment
            in (block_alignments or {}).items()
        },
        "analysis_settings": analysis_settings or {},
        "dotplots": {
            "format_version": 1,
            "pairs": [
                {
                    "pair_id": record.pair_id,
                    "x_sample": record.x_sample,
                    "y_sample": record.y_sample,
                    "svg_rel_path": record.svg_rel_path,
                }
                for record in (dotplots or [])
            ],
        },
    }
