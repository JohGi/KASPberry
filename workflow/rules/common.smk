from pathlib import Path
import sys

sys.path.insert(0, str(Path(workflow.current_basedir) / "../scripts"))

from input_tables import read_genotypes


wildcard_constraints:
    sample="[^/]+",
    pair_id="[^/]+"


# ---------------------------------------------------------------------------
# Project and inputs
# ---------------------------------------------------------------------------

OUTDIR = Path("results")
PROJECT_TITLE = config["project"]["name"]

SCRIPTS_DIR = Path(workflow.current_basedir) / "../scripts"

GENOTYPES_TSV = Path(config["inputs"]["genotypes"])
GENOTYPES = read_genotypes(GENOTYPES_TSV)

GENOTYPE_NAMES = [
    record["genotype"]
    for record in GENOTYPES
]

REGION_FASTA_BY_GENOTYPE = {
    record["genotype"]: record["region_fasta"]
    for record in GENOTYPES
}

NB_GENOTYPES = len(GENOTYPES)


# ---------------------------------------------------------------------------
# KASP genotypes
# ---------------------------------------------------------------------------

KASP_GENOTYPES = [
    record
    for record in GENOTYPES
    if record["genome_fasta"]
]

KASP_GENOTYPE_NAMES = [
    record["genotype"]
    for record in KASP_GENOTYPES
]

KASP_GENOME_FASTA_BY_GENOTYPE = {
    record["genotype"]: record["genome_fasta"]
    for record in KASP_GENOTYPES
}


# ---------------------------------------------------------------------------
# Main output directories
# ---------------------------------------------------------------------------

WORK_DIR = OUTDIR / ".work"

REGIONS_DIR = OUTDIR / "regions"
SNP_DIR = OUTDIR / "snps"
KASP_DIR = OUTDIR / "kasp"
REPORTS_DIR = OUTDIR / "reports"

LOG_DIR = OUTDIR / "logs"
BENCHMARK_DIR = OUTDIR / "benchmarks"


# ---------------------------------------------------------------------------
# Report assets
# ---------------------------------------------------------------------------

REPORT_ASSETS_DIR = REPORTS_DIR / "assets"
DOTPLOT_ASSETS_DIR = REPORT_ASSETS_DIR / "dotplots"
DOTPLOT_GALLERY_ASSETS_DIR = DOTPLOT_ASSETS_DIR / "gallery"
DOTPLOT_VIEWER_ASSETS_DIR = DOTPLOT_ASSETS_DIR / "viewer"

# ---------------------------------------------------------------------------
# Working directories
# ---------------------------------------------------------------------------

CLEAN_FASTA_DIR = WORK_DIR / "clean_fasta"
MULTIFASTA_DIR = CLEAN_FASTA_DIR / "multifasta"

SIBELIAZ_DIR = WORK_DIR / "sibeliaz"

BLOCKS_DIR = WORK_DIR / "blocks"
MASK_CHUNK_DIR = BLOCKS_DIR / "mask_chunks"

BLOCK_FASTA_DIR = WORK_DIR / "block_fastas"
BLOCK_FASTA_SPLIT_DIR = BLOCK_FASTA_DIR / "per_block"

MASKED_DIR = WORK_DIR / "masked_block_fastas"
MASKED_CHUNK_DIR = MASKED_DIR / "chunks"

ALIGN_DIR = WORK_DIR / "alignments"

SNP_WORK_DIR = WORK_DIR / "snps"

DOTPLOT_DIR = WORK_DIR / "dotplots"
DOTPLOT_PAF_DIR = DOTPLOT_DIR / "paf"
DOTPLOT_FORMATTED_DIR = DOTPLOT_DIR / "formatted"
DOTPLOT_PDF_DIR = DOTPLOT_DIR / "pdfs"

RUN_SUMMARY_DIR = WORK_DIR / "run_summary"

VIEWER_DIR = WORK_DIR / "viewer"

MASH_DISTANCES_DIR = WORK_DIR / "mash_distances"

BLOCK_STATS_DIR = WORK_DIR / "block_stats"
UNMASKED_ALIGN_DIR = BLOCK_STATS_DIR / "unmasked_alignments"

KIMURA2P_DISTMAT_DIR = BLOCK_STATS_DIR / "kimura2p_distances"
KIMURA2P_DISTMAT_MATRIX_DIR = KIMURA2P_DISTMAT_DIR / "matrices"
KIMURA2P_DISTMAT_CHUNK_DIR = KIMURA2P_DISTMAT_DIR / "chunks"

KASP_WORK_DIR = WORK_DIR / "kasp"

POLYMARKER_INPUT_DIR = KASP_WORK_DIR / "polymarker_inputs"
POLYMARKER_DIR = KASP_WORK_DIR / "polymarker"
POLYMARKER_SUMMARY_DIR = KASP_WORK_DIR / "polymarker_summary"

IN_SILICO_DIR = KASP_WORK_DIR / "in_silico_validation"
VALIDATION_DIR = KASP_WORK_DIR / "validation"

AGGREGATION_DIR = WORK_DIR / "aggregation"


# ---------------------------------------------------------------------------
# Core FASTA and block files
# ---------------------------------------------------------------------------

CLEAN_FASTAS = expand(
    CLEAN_FASTA_DIR / "{sample}.fasta",
    sample=GENOTYPE_NAMES,
)

ALL_GENOMES_FASTA = MULTIFASTA_DIR / "all_genomes.fasta"

SIBELIAZ_GFF = SIBELIAZ_DIR / "blocks_coords.gff"

FILTERED_GFF = BLOCKS_DIR / "collinear_blocks.gff"
BLOCK_LIST = BLOCKS_DIR / "kept_blocks.list"

ALL_BLOCKS_RAW_FASTA = BLOCK_FASTA_DIR / "all_blocks.raw.fasta"

BLOCK_COORDINATES_TSV = REGIONS_DIR / "collinear_blocks.tsv"


# ---------------------------------------------------------------------------
# SNP files
# ---------------------------------------------------------------------------

SNP_VCF = SNP_DIR / "detected_snps.vcf"
DIAGNOSTIC_SNPS_VCF = SNP_DIR / "diagnostic_snps.vcf"

SNP_POS_LONG_TSV = SNP_DIR / "snp_coordinates.tsv"
SNP_SUMMARY_TSV = SNP_DIR / "snp_summary.tsv"

SNP_POS_WIDE_TSV = SNP_WORK_DIR / "snp_positions_wide.tsv"
DIAGNOSTIC_STATUS_TSV = SNP_WORK_DIR / "diagnostic_status.tsv"

GROUP_A_LIST = SNP_WORK_DIR / "group_a_genotypes.list"
GROUP_B_LIST = SNP_WORK_DIR / "group_b_genotypes.list"


# ---------------------------------------------------------------------------
# Dotplots and viewer files
# ---------------------------------------------------------------------------

DOTPLOT_GALLERY_HTML = REPORTS_DIR / "dotplots.html"

DOTPLOT_MANIFEST = VIEWER_DIR / "dotplots_manifest.json"
GFF_TRACKS_JSON = VIEWER_DIR / "gff_tracks.json"


# ---------------------------------------------------------------------------
# Distance and block statistics
# ---------------------------------------------------------------------------

MASH_MATRIX = MASH_DISTANCES_DIR / "mash.matrix.tsv"

MASKED_BLOCK_N_STATS_TSV = (
    BLOCK_STATS_DIR / "masked_block_n_stats.tsv"
)


# ---------------------------------------------------------------------------
# Run summary
# ---------------------------------------------------------------------------

RUN_SUMMARY_JSON = RUN_SUMMARY_DIR / "run_summary.json"
RUN_SUMMARY_TXT = REPORTS_DIR / "run_summary.txt"


# ---------------------------------------------------------------------------
# PolyMarker intermediate files
# ---------------------------------------------------------------------------

POLYMARKER_DESIGN_STATUS_TSV = (
    POLYMARKER_SUMMARY_DIR / "polymarker_design_status.tsv"
)

POLYMARKER_DESIGN_STATUS_BY_GENOTYPE_TSV = (
    POLYMARKER_SUMMARY_DIR
    / "polymarker_design_status_by_genotype.tsv"
)

POLYMARKER_ASSAYS_TSV = (
    POLYMARKER_SUMMARY_DIR / "polymarker_assays.tsv"
)


# ---------------------------------------------------------------------------
# In silico validation intermediate files
# ---------------------------------------------------------------------------

IN_SILICO_VALIDATION_STATUS_TSV = (
    VALIDATION_DIR / "in_silico_validation_status.tsv"
)

IN_SILICO_ASSAY_STATUS_TSV = (
    VALIDATION_DIR / "in_silico_assay_status.tsv"
)

IN_SILICO_ASSAY_STATUS_BY_GENOTYPE_TSV = (
    VALIDATION_DIR / "in_silico_assay_status_by_genotype.tsv"
)


# ---------------------------------------------------------------------------
# Public KASP outputs
# ---------------------------------------------------------------------------

CANDIDATE_SNPS_VCF = KASP_DIR / "candidate_snps.vcf"

KASP_SUMMARY_TSV = KASP_DIR / "kasp_snp_summary.tsv"
ASSAY_SUMMARY_TSV = KASP_DIR / "assay_summary.tsv"
CANDIDATE_ASSAYS_TSV = KASP_DIR / "candidate_assays.tsv"
PRIMERS_TO_ORDER_TSV = KASP_DIR / "primers_to_order.tsv"


# ---------------------------------------------------------------------------
# Repeat masking configuration
# ---------------------------------------------------------------------------

te_lib_value = config["snps"]["repeat_masking"]["library"]

TE_LIB = Path(te_lib_value) if te_lib_value else None
ALIGNMENT_FASTA_DIR = MASKED_DIR


# ---------------------------------------------------------------------------
# Dynamic checkpoint helpers
# ---------------------------------------------------------------------------

def get_split_block_dir(_wildcards=None) -> Path:
    """Return the checkpoint output directory containing per-block FASTA files."""
    return Path(
        checkpoints.split_block_fastas.get().output[0]
    )


def get_chunk_ids(_wildcards=None) -> list[str]:
    """Return all chunk IDs after checkpoint completion."""
    chunk_dir = Path(
        checkpoints.split_block_list_into_chunks.get().output.chunk_dir
    )

    return sorted(
        path.stem
        for path in chunk_dir.glob("*.list")
    )


def get_alignment_inputs(wildcards) -> list[Path]:
    """Return prerequisite inputs for one alignment chunk."""
    return [
        MASK_CHUNK_DIR / f"{wildcards.chunk_id}.list",
        MASKED_CHUNK_DIR / f"{wildcards.chunk_id}.done",
    ]


def get_align_chunk_sentinels(_wildcards=None) -> list[Path]:
    """Return all alignment chunk completion markers."""
    return [
        ALIGN_DIR / f"{chunk_id}.done"
        for chunk_id in get_chunk_ids()
    ]


def get_unmasked_align_chunk_sentinels(
    _wildcards=None,
) -> list[Path]:
    """Return all unmasked alignment chunk completion markers."""
    return [
        UNMASKED_ALIGN_DIR / f"{chunk_id}.done"
        for chunk_id in get_chunk_ids()
    ]


def get_distmat_chunk_sentinels(
    _wildcards=None,
) -> list[Path]:
    """Return all distance-matrix chunk completion markers."""
    return [
        KIMURA2P_DISTMAT_CHUNK_DIR / f"{chunk_id}.done"
        for chunk_id in get_chunk_ids()
    ]
