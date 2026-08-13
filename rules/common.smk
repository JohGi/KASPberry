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
PROJECT_TITLE = config.get("project", {}).get("name", "Project")

CLEAN_FASTA_DIR = OUTDIR / "01_clean_fasta"
MULTIFASTA_DIR = CLEAN_FASTA_DIR / "multifasta"
SIBELIAZ_DIR = OUTDIR / "02_sibeliaz"
FILTERED_BLOCKS_DIR = OUTDIR / "03_filtered_blocks"
MASK_CHUNK_DIR = FILTERED_BLOCKS_DIR / "mask_chunks"
BLOCK_FASTA_DIR = OUTDIR / "04_block_fastas"
BLOCK_FASTA_SPLIT_DIR = BLOCK_FASTA_DIR / "per_block"
MASKED_DIR = OUTDIR / "05_masked_block_fastas"
MASKED_CHUNK_DIR = MASKED_DIR / "chunks"
ALIGN_DIR = OUTDIR / "06_alignments"
SNP_DIR = OUTDIR / "07_snps"
FILTERED_SNP_DIR = OUTDIR / "08_filtered_snps"
SNP_POS_DIR = OUTDIR / "09_snp_positions"
DOTPLOT_DIR = OUTDIR / "10_dotplots"
DOTPLOT_PAF_DIR = DOTPLOT_DIR / "paf"
DOTPLOT_FORMATTED_DIR = DOTPLOT_DIR / "formatted"
DOTPLOT_PDF_DIR = DOTPLOT_DIR / "pdfs"
DOTPLOT_SVG_DIR = DOTPLOT_DIR / "svgs"
DOTPLOT_COMBINED_DIR = DOTPLOT_DIR / "combined"
SUMMARY_STATS_DIR = OUTDIR / "11_summary_stats"
REGION_TRACK_DIR = OUTDIR / "12_region_viewer"
DOTPLOT_ONLY_SVG_DIR = REGION_TRACK_DIR / "viewer_assets/dotplots"
MASH_DISTANCES_DIR = OUTDIR / "13_mash_distances"
BLOCK_STATS_DIR = OUTDIR / "14_block_stats"
UNMASKED_ALIGN_DIR = BLOCK_STATS_DIR / "unmasked_alignments"
KIMURA2P_DISTMAT_DIR = BLOCK_STATS_DIR / "kimura2p_distances"
KIMURA2P_DISTMAT_MATRIX_DIR = KIMURA2P_DISTMAT_DIR / "matrices"
KIMURA2P_DISTMAT_CHUNK_DIR = KIMURA2P_DISTMAT_DIR / "chunks"
LOG_DIR = OUTDIR / "logs"
BENCHMARK_DIR = OUTDIR / "benchmarks"

CLEAN_FASTAS = expand(CLEAN_FASTA_DIR / "{sample}.fasta", sample=GENOTYPE_NAMES)
ALL_GENOMES_FASTA = MULTIFASTA_DIR / "all_genomes.fasta"
SIBELIAZ_GFF = SIBELIAZ_DIR / "blocks_coords.gff"
FILTERED_GFF = FILTERED_BLOCKS_DIR / "filtered_blocks.gff"
BLOCK_LIST = FILTERED_BLOCKS_DIR / "kept_blocks.list"
BLOCK_COORDINATES_TSV = FILTERED_BLOCKS_DIR / "block_coords.tsv"
ALL_BLOCKS_RAW_FASTA = BLOCK_FASTA_DIR / "all_blocks.raw.fasta"
SNP_VCF = SNP_DIR / "snps.vcf"
GROUP_A_LIST = FILTERED_SNP_DIR / "group_a_samples.list"
GROUP_B_LIST = FILTERED_SNP_DIR / "group_b_samples.list"
FILTERED_SNP_VCF = FILTERED_SNP_DIR / "filtered_snps.vcf"
SNP_POS_LONG_TSV = SNP_POS_DIR / "snp_positions_long.tsv"
SNP_POS_WIDE_TSV = SNP_POS_DIR / "snp_positions_wide.tsv"
DOTPLOT_GALLERY_HTML = DOTPLOT_COMBINED_DIR / "dotplots_gallery.html"
SUMMARY_STATS_JSON = SUMMARY_STATS_DIR / "summary_stats.json"
SUMMARY_STATS_TXT = SUMMARY_STATS_DIR / "summary_stats.txt"
DOTPLOT_MANIFEST = REGION_TRACK_DIR / "dotplots_manifest.json"
REGION_TRACK_HTML = REGION_TRACK_DIR / "region_tracks.html"
MASHTREE_MATRIX = MASH_DISTANCES_DIR / "mashtree.matrix.tsv"
MASHTREE_TREE = MASH_DISTANCES_DIR / "mashtree.dnd"
MASKED_BLOCK_N_STATS_TSV = BLOCK_STATS_DIR / "masked_block_n_stats.tsv"
GFF_TRACKS_JSON = REGION_TRACK_DIR / "gff_tracks.json"

NB_GENOTYPES = len(GENOTYPES)
te_lib_value = config.get("snps", {}).get("repeat_masking", {}).get("library", "")
TE_LIB = Path(te_lib_value) if te_lib_value else None
USE_MASKING = TE_LIB is not None

ALIGNMENT_FASTA_DIR = MASKED_DIR if USE_MASKING else BLOCK_FASTA_SPLIT_DIR
