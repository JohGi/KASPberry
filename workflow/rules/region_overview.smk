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
        REPORTS_DIR / f"region_viewer_{mode}.html"
    ]


def get_viewer_snp_summary(wildcards):
    """Return the aggregated SNP table for the requested viewer mode."""
    return {
        "snps": SNP_SUMMARY_TSV,
        "kasp": KASP_SUMMARY_TSV,
    }[wildcards.mode]


def get_viewer_assay_summary(wildcards):
    """Return the assay summary for KASP mode only."""
    if wildcards.mode == "kasp":
        return ASSAY_SUMMARY_TSV

    return []


def get_analysis_settings(mode):
    """Return analysis settings effectively used by the requested workflow mode."""
    settings = {
        "minimum_block_length_bp": config["snps"]["min_block_length"],
        "snp_groups": {
            group: [
                record["genotype"]
                for record in GENOTYPES
                if record["group"] == group
            ]
            for group in sorted(
                {
                    record["group"]
                    for record in GENOTYPES
                    if record["group"]
                }
            )
        },
        "repeat_masking": {
            "simple_repeats_and_low_complexity": True,
            "custom_repeat_library": config["snps"]["repeat_masking"][
                "library"
            ],
        },
    }

    advanced_options = {
        "sibeliaz": config["advanced"]["sibeliaz"]["extra_options"],
        "mafft": config["advanced"]["alignment"]["mafft_options"],
    }

    if mode == "kasp":
        settings["kasp_assay_design"] = {
            "polymarker_subgenomes": config["kasp"][
                "polymarker_genomes"
            ],
            "genotypes": KASP_GENOTYPE_NAMES,
            "mfeprimer_min_tm": config["kasp"]["mfeprimer"][
                "specificity"
            ]["min_tm"],
            "mfeprimer_dimer_max_dg": config["kasp"]["mfeprimer"][
                "dimer"
            ]["max_dg"],
        }

        advanced_options.update(
            {
                "mfeprimer_specificity": config["advanced"][
                    "mfeprimer"
                ]["specificity_extra_options"],
                "mfeprimer_dimer": config["advanced"]["mfeprimer"][
                    "dimer_extra_options"
                ],
                "mfeprimer_hairpin": config["advanced"]["mfeprimer"][
                    "hairpin_extra_options"
                ],
            }
        )

    non_empty_advanced_options = {
        name: options
        for name, options in advanced_options.items()
        if options
    }

    if non_empty_advanced_options:
        settings["advanced_options"] = non_empty_advanced_options

    return settings


rule write_analysis_settings_json:
    wildcard_constraints:
        mode="snps|kasp"

    output:
        settings=ANALYSIS_SETTINGS_JSON

    params:
        settings=lambda wildcards: get_analysis_settings(wildcards.mode)

    benchmark:
        BENCHMARK_DIR / "write_analysis_settings_json/{mode}.tsv"

    log:
        stdout=LOG_DIR / "write_analysis_settings_json/{mode}.stdout",
        stderr=LOG_DIR / "write_analysis_settings_json/{mode}.stderr"

    run:
        import json

        output_path = Path(output.settings)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        with output_path.open("w", encoding="utf-8") as handle:
            json.dump(params.settings, handle, indent=2)
            handle.write("\n")


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
    conda:
        "../envs/workflow-runtime.yaml"
    input:
        svgs=DOTPLOT_VIEWER_SVGS
    output:
        manifest=DOTPLOT_MANIFEST
    benchmark:
        BENCHMARK_DIR / "build_dotplot_manifest" / "dotplots_manifest.tsv"
    log:
        stdout=LOG_DIR / "build_dotplot_manifest" / "dotplots_manifest.stdout",
        stderr=LOG_DIR / "build_dotplot_manifest" / "dotplots_manifest.stderr"
    shell:
        r"""
        mkdir -p "$(dirname "{log.stdout}")"

        python3 "{SCRIPTS_DIR}/build_dotplot_manifest.py" \
            --genotypes "{GENOTYPES_TSV}" \
            --svg-dir "{DOTPLOT_VIEWER_ASSETS_DIR}" \
            --output "{output.manifest}" \
            1> "{log.stdout}" \
            2> "{log.stderr}"
        """

rule generate_region_viewer:
    conda:
        "../envs/workflow-runtime.yaml"
    wildcard_constraints:
        mode="snps|kasp"

    input:
        samples_tsv=GENOTYPES_TSV,
        block_coords_tsv=BLOCK_COORDINATES_TSV,
        snp_long=SNP_POS_LONG_TSV,
        snp_summary=get_viewer_snp_summary,
        assay_summary=get_viewer_assay_summary,
        fastas=CLEAN_FASTAS,
        stats_json=RUN_SUMMARY_JSON,
        mash_dists_tsv=MASH_MATRIX,
        n_stats_tsv=MASKED_BLOCK_N_STATS_TSV,
        align_sentinels=get_align_chunk_sentinels,
        distmat_sentinels=get_distmat_chunk_sentinels,
        gff_tracks_json=GFF_TRACKS_JSON,
        dotplot_manifest=DOTPLOT_MANIFEST,
        analysis_settings_json=ANALYSIS_SETTINGS_JSON,

    output:
        html=REPORTS_DIR / "region_viewer_{mode}.html"

    params:
        mode=lambda wildcards: wildcards.mode,

    benchmark:
        BENCHMARK_DIR / "generate_region_viewer/{mode}.tsv"

    log:
        stdout=LOG_DIR / "generate_region_viewer/{mode}.stdout",
        stderr=LOG_DIR / "generate_region_viewer/{mode}.stderr"

    shell:
        r"""
        mkdir -p "$(dirname "{output.html}")" "$(dirname "{log.stdout}")"

        assay_args=()

        if [[ "{params.mode}" == "kasp" ]]; then
            assay_args=(
                --assay-summary "{input.assay_summary}"
            )
        fi

        python3 "{SCRIPTS_DIR}/generate_region_viewer.py" \
            --mode "{params.mode}" \
            --samples-tsv "{input.samples_tsv}" \
            --block-coords-tsv "{input.block_coords_tsv}" \
            --snp-long "{input.snp_long}" \
            --snp-summary "{input.snp_summary}" \
            "${{assay_args[@]}}" \
            --fasta-dir "{CLEAN_FASTA_DIR}" \
            --summary-stats-json "{input.stats_json}" \
            --mash-matrix "{input.mash_dists_tsv}" \
            --kimura2p-distmat-dir "{KIMURA2P_DISTMAT_MATRIX_DIR}" \
            --masked-align-dir "{ALIGN_DIR}" \
            --masked-block-n-stats "{input.n_stats_tsv}" \
            --gff-tracks-json "{input.gff_tracks_json}" \
            --dotplot-manifest-json {input.dotplot_manifest} \
            --analysis-settings-json "{input.analysis_settings_json}" \
            --title "{PROJECT_TITLE}" \
            --output "{output.html}" \
            1> "{log.stdout}" \
            2> "{log.stderr}"
        """
