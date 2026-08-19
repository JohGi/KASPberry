from pathlib import Path
from input_tables import read_annotations


def get_gff_tracks(config):
    """Return configured GFF tracks."""
    annotation_path = config["inputs"]["annotations"]

    if not annotation_path:
        return {}

    return read_annotations(Path(annotation_path))


def get_gff_track_files(gff_tracks):
    """Return all configured GFF track files."""
    return [
        gff_path
        for sample_tracks in gff_tracks.values()
        for gff_path in sample_tracks.values()
    ]


GFF_TRACKS = get_gff_tracks(config)
GFF_TRACK_FILES = get_gff_track_files(GFF_TRACKS)


def get_region_viewer_outputs(mode):
    """Return the region viewer HTML output for one workflow mode."""
    return [
        REGION_TRACK_DIR / f"region_tracks.{mode}.html"
    ]


def get_viewer_snp_summary(wildcards):
    """Return the aggregated SNP table for the requested viewer mode."""
    return {
        "snps": SNP_SUMMARY_TSV,
        "kasp": KASP_SUMMARY_TSV,
    }[wildcards.mode]


rule write_gff_tracks_json:
    input:
        gff_files=GFF_TRACK_FILES
    output:
        GFF_TRACKS_JSON
    benchmark:
        BENCHMARK_DIR / "write_gff_tracks_json.tsv"
    log:
        stdout=LOG_DIR / "write_gff_tracks_json" / "write_gff_tracks_json.stdout",
        stderr=LOG_DIR / "write_gff_tracks_json" / "write_gff_tracks_json.stderr"
    run:
        import json

        output_path = Path(output[0])
        output_path.parent.mkdir(parents=True, exist_ok=True)

        with output_path.open("w", encoding="utf-8") as handle:
            json.dump(GFF_TRACKS, handle, indent=2)
            handle.write("\n")

rule build_dotplot_manifest:
    input:
        svgs=DOTPLOT_ONLY_SVGS
    output:
        manifest=DOTPLOT_MANIFEST
    benchmark:
        BENCHMARK_DIR / "build_dotplot_manifest" / "dotplots_manifest.tsv"
    log:
        stdout=LOG_DIR / "build_dotplot_manifest" / "dotplots_manifest.stdout",
        stderr=LOG_DIR / "build_dotplot_manifest" / "dotplots_manifest.stderr"
    shell:
        r"""
        mkdir -p "{DOTPLOT_COMBINED_DIR}" "$(dirname "{log.stdout}")"

        python3 "{SCRIPTS_DIR}/build_dotplot_manifest.py" \
            --genotypes "{GENOTYPES_TSV}" \
            --svg-dir "{DOTPLOT_ONLY_SVG_DIR}" \
            --output "{output.manifest}" \
            1> "{log.stdout}" \
            2> "{log.stderr}"
        """

rule generate_region_viewer:
    wildcard_constraints:
        mode="snps|kasp"

    input:
        samples_tsv=GENOTYPES_TSV,
        block_coords_tsv=BLOCK_COORDINATES_TSV,
        snp_long=SNP_POS_LONG_TSV,
        snp_summary=get_viewer_snp_summary,
        fastas=CLEAN_FASTAS,
        stats_json=SUMMARY_STATS_JSON,
        mash_dists_tsv=MASHTREE_MATRIX,
        n_stats_tsv=MASKED_BLOCK_N_STATS_TSV,
        align_sentinels=get_align_chunk_sentinels,
        distmat_sentinels=get_distmat_chunk_sentinels,
        gff_tracks_json=GFF_TRACKS_JSON,
        dotplot_manifest=DOTPLOT_MANIFEST,

    output:
        html=REGION_TRACK_DIR / "region_tracks.{mode}.html"

    benchmark:
        BENCHMARK_DIR / "generate_region_viewer/{mode}.tsv"

    log:
        stdout=LOG_DIR / "generate_region_viewer/{mode}.stdout",
        stderr=LOG_DIR / "generate_region_viewer/{mode}.stderr"

    shell:
        r"""
        mkdir -p "{REGION_TRACK_DIR}" "$(dirname "{log.stdout}")"

        python3 "{SCRIPTS_DIR}/generate_region_viewer.py" \
            --samples-tsv "{input.samples_tsv}" \
            --block-coords-tsv "{input.block_coords_tsv}" \
            --snp-long "{input.snp_long}" \
            --fasta-dir "{CLEAN_FASTA_DIR}" \
            --summary-stats-json "{input.stats_json}" \
            --mash-matrix "{input.mash_dists_tsv}" \
            --kimura2p-distmat-dir "{KIMURA2P_DISTMAT_MATRIX_DIR}" \
            --masked-align-dir "{ALIGN_DIR}" \
            --masked-block-n-stats "{input.n_stats_tsv}" \
            --gff-tracks-json "{input.gff_tracks_json}" \
            --dotplot-manifest-json {input.dotplot_manifest} \
            --config-yaml "{workflow.configfiles[0]}" \
            --title "{PROJECT_TITLE}" \
            --output "{output.html}" \
            1> "{log.stdout}" \
            2> "{log.stderr}"
        """
