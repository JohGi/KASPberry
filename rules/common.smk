from pathlib import Path
import sys
sys.path.insert(0, str(Path(workflow.current_basedir) / "../scripts"))
from input_tables import read_genotypes

wildcard_constraints:
    sample="[^/]+",
    pair_id="[^/]+"


def get_split_block_dir(_wildcards=None) -> Path:
    """Return the checkpoint output directory containing per-block FASTA files."""
    return Path(checkpoints.split_block_fastas.get().output[0])


def get_chunk_ids(_wildcards=None) -> list[str]:
    """Return all chunk IDs after checkpoint completion."""
    chunk_dir = Path(checkpoints.split_block_list_into_chunks.get().output.chunk_dir)
    return sorted(path.stem for path in chunk_dir.glob("*.list"))


def get_alignment_inputs(wildcards) -> list[Path]:
    """Return prerequisite inputs for one alignment chunk."""
    inputs = [MASK_CHUNK_DIR / f"{wildcards.chunk_id}.list"]

    if USE_MASKING:
        inputs.append(MASKED_CHUNK_DIR / f"{wildcards.chunk_id}.done")
    else:
        inputs.append(get_split_block_dir())

    return inputs


def get_align_chunk_sentinels(_wildcards=None) -> list[Path]:
    """Return all alignment chunk completion markers after checkpoint completion."""
    return [ALIGN_DIR / f"{chunk_id}.done" for chunk_id in get_chunk_ids()]


def get_unmasked_align_chunk_sentinels(_wildcards=None) -> list[Path]:
    """Return all unmasked alignment chunk completion markers after checkpoint completion."""
    return [UNMASKED_ALIGN_DIR / f"{chunk_id}.done" for chunk_id in get_chunk_ids()]


def get_distmat_chunk_sentinels(_wildcards=None) -> list[Path]:
    """Return all distmat chunk completion markers after checkpoint completion."""
    return [
        KIMURA2P_DISTMAT_CHUNK_DIR / f"{chunk_id}.done"
        for chunk_id in get_chunk_ids()
    ]


GENOTYPES_TSV = Path(config["inputs"]["genotypes"])
GENOTYPES = read_genotypes(GENOTYPES_TSV)
GENOTYPE_NAMES = [record["genotype"] for record in GENOTYPES]
REGION_FASTA_BY_GENOTYPE = {
    record["genotype"]: record["region_fasta"] for record in GENOTYPES
}

OUTDIR = Path(config["project"]["output_dir"])
SCRIPTS_DIR = Path(workflow.current_basedir) / "../scripts"
PROJECT_TITLE = config["project"]["name"]

CLEAN_FASTA_DIR = OUTDIR / ".work/clean_fasta"
MULTIFASTA_DIR = CLEAN_FASTA_DIR / "multifasta"
SIBELIAZ_DIR = OUTDIR / ".work/sibeliaz"
FILTERED_BLOCKS_DIR = OUTDIR / ".work/blocks"
MASK_CHUNK_DIR = FILTERED_BLOCKS_DIR / "mask_chunks"
BLOCK_FASTA_DIR = OUTDIR / ".work/block_fastas"
BLOCK_FASTA_SPLIT_DIR = BLOCK_FASTA_DIR / "per_block"
MASKED_DIR = OUTDIR / ".work/masked_block_fastas"
MASKED_CHUNK_DIR = MASKED_DIR / "chunks"
ALIGN_DIR = OUTDIR / ".work/alignments"
SNP_DIR = OUTDIR / "snps"
FILTERED_SNP_DIR = OUTDIR / ".work/snps"
SNP_POS_DIR = OUTDIR / ".work/snps"
DOTPLOT_DIR = OUTDIR / ".work/dotplots"
DOTPLOT_PAF_DIR = DOTPLOT_DIR / "paf"
DOTPLOT_FORMATTED_DIR = DOTPLOT_DIR / "formatted"
DOTPLOT_PDF_DIR = DOTPLOT_DIR / "pdfs"
DOTPLOT_SVG_DIR = DOTPLOT_DIR / "svgs"
DOTPLOT_COMBINED_DIR = DOTPLOT_DIR / "combined"
RUN_SUMMARY_DIR = OUTDIR / ".work/run_summary"
REGION_TRACK_DIR = OUTDIR / ".work/viewer"
DOTPLOT_ONLY_SVG_DIR = REGION_TRACK_DIR / "viewer_assets/dotplots"
MASH_DISTANCES_DIR = OUTDIR / ".work/mash_distances"
BLOCK_STATS_DIR = OUTDIR / ".work/block_stats"
UNMASKED_ALIGN_DIR = BLOCK_STATS_DIR / "unmasked_alignments"
KIMURA2P_DISTMAT_DIR = BLOCK_STATS_DIR / "kimura2p_distances"
KIMURA2P_DISTMAT_MATRIX_DIR = KIMURA2P_DISTMAT_DIR / "matrices"
KIMURA2P_DISTMAT_CHUNK_DIR = KIMURA2P_DISTMAT_DIR / "chunks"
POLYMARKER_INPUT_DIR = OUTDIR / ".work/kasp/polymarker_inputs"
POLYMARKER_DIR = OUTDIR / ".work/kasp/polymarker"
POLYMARKER_SUMMARY_DIR = OUTDIR / ".work/kasp/polymarker_summary"
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
IN_SILICO_DIR = OUTDIR / ".work/kasp/in_silico_validation"
VALIDATION_DIR = OUTDIR / ".work/kasp/validation"
AGGREGATION_DIR = OUTDIR / ".work/aggregation"
LOG_DIR = OUTDIR / "logs"
BENCHMARK_DIR = OUTDIR / "benchmarks"

CLEAN_FASTAS = expand(CLEAN_FASTA_DIR / "{sample}.fasta", sample=GENOTYPE_NAMES)
ALL_GENOMES_FASTA = MULTIFASTA_DIR / "all_genomes.fasta"
SIBELIAZ_GFF = SIBELIAZ_DIR / "blocks_coords.gff"
FILTERED_GFF = OUTDIR / ".work/blocks/collinear_blocks.gff"
BLOCK_LIST = FILTERED_BLOCKS_DIR / "kept_blocks.list"
BLOCK_COORDINATES_TSV = OUTDIR / "regions/collinear_blocks.tsv"
ALL_BLOCKS_RAW_FASTA = BLOCK_FASTA_DIR / "all_blocks.raw.fasta"
SNP_VCF = SNP_DIR / "detected_snps.vcf"
GROUP_A_LIST = FILTERED_SNP_DIR / "group_a_samples.list"
GROUP_B_LIST = FILTERED_SNP_DIR / "group_b_samples.list"
SNP_POS_LONG_TSV = OUTDIR / "snps/snp_coordinates.tsv"
SNP_POS_WIDE_TSV = SNP_POS_DIR / "snp_positions_wide.tsv"
DOTPLOT_GALLERY_HTML = OUTDIR / "reports/dotplots.html"
RUN_SUMMARY_JSON = RUN_SUMMARY_DIR / "run_summary.json"
RUN_SUMMARY_TXT = OUTDIR / "reports/run_summary.txt"
DOTPLOT_MANIFEST = REGION_TRACK_DIR / "dotplots_manifest.json"
MASHTREE_MATRIX = MASH_DISTANCES_DIR / "mashtree.matrix.tsv"
MASHTREE_TREE = MASH_DISTANCES_DIR / "mashtree.dnd"
MASKED_BLOCK_N_STATS_TSV = BLOCK_STATS_DIR / "masked_block_n_stats.tsv"
GFF_TRACKS_JSON = REGION_TRACK_DIR / "gff_tracks.json"

DIAGNOSTIC_STATUS_TSV = FILTERED_SNP_DIR / "diagnostic_status.tsv"
POLYMARKER_DESIGN_STATUS_TSV = POLYMARKER_SUMMARY_DIR / "polymarker_design_status.tsv"
IN_SILICO_VALIDATION_STATUS_TSV = VALIDATION_DIR / "in_silico_validation_status.tsv"
POLYMARKER_ASSAYS_TSV = POLYMARKER_SUMMARY_DIR / "polymarker_assays.tsv"
POLYMARKER_DESIGN_STATUS_BY_GENOTYPE_TSV = POLYMARKER_SUMMARY_DIR / "polymarker_design_status_by_genotype.tsv"
IN_SILICO_ASSAY_STATUS_TSV = VALIDATION_DIR / "in_silico_assay_status.tsv"
IN_SILICO_ASSAY_STATUS_BY_GENOTYPE_TSV = VALIDATION_DIR / "in_silico_assay_status_by_genotype.tsv"
DIAGNOSTIC_SNPS_VCF = OUTDIR / "snps/diagnostic_snps.vcf"
CANDIDATE_SNPS_VCF = OUTDIR / "kasp/candidate_snps.vcf"
SNP_SUMMARY_TSV = OUTDIR / "snps/snp_summary.tsv"
KASP_SUMMARY_TSV = OUTDIR / "kasp/kasp_snp_summary.tsv"
ASSAY_SUMMARY_TSV = OUTDIR / "kasp/assay_summary.tsv"
CANDIDATE_ASSAYS_TSV = OUTDIR / "kasp/candidate_assays.tsv"
PRIMERS_TO_ORDER_TSV = OUTDIR / "kasp/primers_to_order.tsv"

NB_GENOTYPES = len(GENOTYPES)
te_lib_value = config["snps"]["repeat_masking"]["library"]
TE_LIB = Path(te_lib_value) if te_lib_value else None
USE_MASKING = TE_LIB is not None

ALIGNMENT_FASTA_DIR = MASKED_DIR if USE_MASKING else BLOCK_FASTA_SPLIT_DIR
