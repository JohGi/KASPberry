from snakemake.utils import validate

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
include: "rules/in_silico_validation.smk"

SNP_TARGET_OUTPUTS = [
    SNP_POS_LONG_TSV,
    SNP_POS_WIDE_TSV,
    *get_region_viewer_outputs(),
    DOTPLOT_GALLERY_HTML,
    SUMMARY_STATS_TXT,
]

KASP_TARGET_OUTPUTS = [
    VALIDATION_DIR / "assay_status.tsv",
    VALIDATION_DIR / "assay_status_by_genotype.tsv",
    VALIDATION_DIR / "snp_status.tsv",
    VALIDATION_DIR / "validated_assays.tsv",
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
