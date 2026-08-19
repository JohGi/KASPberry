from snakemake.utils import validate

validate(config, "workflow/schemas/config.schema.yaml")

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
include: "rules/in_silico_validation.smk"
include: "rules/aggregation.smk"

SNP_TARGET_OUTPUTS = [
    SNP_POS_LONG_TSV,
    SNP_POS_WIDE_TSV,
    SNP_SUMMARY_TSV,
    DIAGNOSTIC_SNPS_VCF,
    *get_region_viewer_outputs("snps"),
    DOTPLOT_GALLERY_HTML,
    SUMMARY_STATS_TXT,
]

KASP_TARGET_OUTPUTS = [
    SNP_POS_LONG_TSV,
    SNP_POS_WIDE_TSV,
    KASP_SUMMARY_TSV,
    VALIDATED_SNPS_VCF,
    ASSAY_SUMMARY_TSV,
    VALIDATED_ASSAYS_TSV,
    PRIMERS_TO_ORDER_TSV,
    *get_region_viewer_outputs("kasp"),
    DOTPLOT_GALLERY_HTML,
    SUMMARY_STATS_TXT,
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
