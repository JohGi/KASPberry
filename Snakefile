from snakemake.utils import validate

# configfile: "config/config.yaml"

validate(config, "workflow/schemas/config.schema.yaml", set_default=False)

include: "rules/common.smk"
include: "rules/mash_dists.smk"
include: "rules/blocks.smk"
include: "rules/snps.smk"
include: "rules/coordinate_mapping.smk"
include: "rules/dotplots.smk"
include: "rules/summary_stats.smk"
include: "rules/block_stats.smk"
include: "rules/region_overview.smk"
include: "rules/kasp.smk"


SNP_TARGET_OUTPUTS = [
    SNP_POS_LONG_TSV,
    SNP_POS_WIDE_TSV,
    *get_region_viewer_outputs(),
    DOTPLOT_GALLERY_HTML,
    SUMMARY_STATS_TXT,
]

KASP_TARGET_OUTPUTS = [
    *SNP_TARGET_OUTPUTS,
    POLYMARKER_SUMMARY_DIR / "polymarker_snp_status.tsv",
    POLYMARKER_SUMMARY_DIR / "polymarker_snp_status_long.tsv",
    POLYMARKER_SUMMARY_DIR / "polymarker_assays.tsv",
    POLYMARKER_SUMMARY_DIR / "primers.tsv",
    POLYMARKER_SUMMARY_DIR / "primers_control_pairs.tsv",
    POLYMARKER_SUMMARY_DIR / "primers_with_tails.fasta",
]


rule all:
    input:
        SNP_TARGET_OUTPUTS

rule snps:
    input:
        SNP_TARGET_OUTPUTS

rule kasp:
    input:
        KASP_TARGET_OUTPUTS
