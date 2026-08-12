from snakemake.utils import validate
from pathlib import Path
import csv
import re
import polars as pl
from itertools import combinations
import sys
sys.path.insert(0, str(Path(workflow.current_basedir) / "../scripts"))
from genotypes import read_genotypes, read_annotations
from input_validation import (
    check_genotypes_table,
    check_annotations_table,
)

wildcard_constraints:
    sample="[^/]+",
    pair_id="[^/]+"


def resolve_snp_filter_groups(
    genotypes: list[dict[str, str]],
) -> tuple[list[str], list[str], bool]:
    """Resolve SNP filtering groups from the genotype table."""
    groups = sorted({row["group"] for row in genotypes if row.get("group")})
    if len(groups) != 2:
        return [], [], False

    group_a = [row["genotype"] for row in genotypes if row.get("group") == groups[0]]
    group_b = [row["genotype"] for row in genotypes if row.get("group") == groups[1]]

    return group_a, group_b, True

def resolve_dotplot_pairs(
    sample_names: list[str],
    config: dict,
) -> list[tuple[str, str]]:
    """Resolve pairwise dotplot comparisons from config."""
    pivot = str(config.get("viewer", {}).get("dotplot_reference") or "").strip()

    if not pivot:
        return list(combinations(sample_names, 2))

    if pivot not in sample_names:
        raise ValueError(
            f"Unknown viewer.dotplot_reference: {pivot!r}. "
            f"Expected one of: {sample_names}"
        )

    return [(pivot, sample) for sample in sample_names if sample != pivot]


def build_pair_id(sample_a: str, sample_b: str) -> str:
    """Build a stable pair identifier."""
    return f"{sample_a}__vs__{sample_b}"


def split_pair_id(pair_id: str) -> tuple[str, str]:
    """Decode a pair identifier."""
    parts = pair_id.split("__vs__")
    if len(parts) != 2 or not parts[0] or not parts[1]:
        raise ValueError(f"Invalid pair_id: {pair_id!r}")
    return parts[0], parts[1]


def get_pair_sample_a(wildcards) -> str:
    """Return sample A for a pair wildcard."""
    sample_a, _sample_b = split_pair_id(wildcards.pair_id)
    return sample_a


def get_pair_sample_b(wildcards) -> str:
    """Return sample B for a pair wildcard."""
    _sample_a, sample_b = split_pair_id(wildcards.pair_id)
    return sample_b


def get_gff_tracks(config, genotype_names):
    """Return configured GFF tracks after validating their structure."""

    annotation_path = config.get("inputs", {}).get("annotations")

    if not annotation_path:
        return {}

    annotation_path = Path(annotation_path)

    annotations_df = pl.read_csv(
        annotation_path,
        separator="\t",
        null_values="",
    )

    validate(
        annotations_df,
        SCHEMA_DIR / "annotations.schema.yaml",
        set_default=False,
    )

    check_annotations_table(
        annotations_df,
        genotype_names=set(genotype_names),
    )

    return read_annotations(annotation_path, set(genotype_names))


def get_gff_track_files(gff_tracks):
    """Return all configured GFF track files."""
    return [
        gff_path
        for sample_tracks in gff_tracks.values()
        for gff_path in sample_tracks.values()
    ]


def get_split_block_dir(_wildcards=None) -> Path:
    """Return the checkpoint output directory containing per-block FASTA files."""
    return Path(checkpoints.split_block_fastas.get().output[0])


def get_chunk_list_dir(_wildcards=None) -> Path:
    """Return the checkpoint output directory containing chunk list files."""
    return Path(checkpoints.split_block_list_into_chunks.get().output.chunk_dir)


def get_chunk_ids(_wildcards=None) -> list[str]:
    """Return all chunk IDs after checkpoint completion."""
    chunk_dir = get_chunk_list_dir()
    return sorted(path.stem for path in chunk_dir.glob("*.list"))

def get_masked_chunk_done_outputs(_wildcards=None) -> list[Path]:
    """Return all masking chunk completion markers after checkpoint completion."""
    return [MASKED_CHUNK_DIR / f"{chunk_id}.done" for chunk_id in get_chunk_ids()]


def get_alignment_fasta_dir() -> str:
    """Return the FASTA directory path used for alignment."""
    if USE_MASKING:
        return str(MASKED_DIR)
    return str(BLOCK_FASTA_SPLIT_DIR)


def get_alignment_inputs(wildcards) -> list[Path]:
    """Return prerequisite inputs for one alignment chunk."""
    inputs = [MASK_CHUNK_DIR / f"{wildcards.chunk_id}.list"]

    if USE_MASKING:
        inputs.append(MASKED_CHUNK_DIR / f"{wildcards.chunk_id}.done")
    else:
        inputs.append(get_split_block_dir())

    return inputs


def get_unmasked_chunk_list(wildcards) -> Path:
    """Return the chunk list for one unmasked alignment chunk."""
    return MASK_CHUNK_DIR / f"{wildcards.chunk_id}.list"


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


def get_final_snp_output() -> Path:
    """Return the final SNP output path depending on filtering settings."""
    if USE_SNP_GROUP_FILTERING:
        return FILTERED_SNP_VCF
    return SNP_VCF


def write_lines(output_path: Path, values: list[str]) -> None:
    """Write one value per line to a text file."""
    with open(output_path, "w", encoding="utf-8") as handle:
        for value in values:
            handle.write(f"{value}\n")


def slugify_marker_set_name(name: str) -> str:
    """Convert a marker set name into a safe filename fragment."""
    slug = re.sub(r"[^A-Za-z0-9._-]+", "_", name.strip())
    slug = re.sub(r"_+", "_", slug)
    return slug.strip("_")


def get_selected_region_viewer_title(wildcards) -> str:
    """Return the title for one selected marker region viewer."""
    marker_set = SELECTED_MARKER_SETS_BY_SLUG[wildcards.marker_set]
    return marker_set["title"]


def get_region_viewer_outputs():
    """Return all region viewer HTML outputs."""
    return [REGION_TRACK_HTML, *REGION_TRACK_SELECTED_HTMLS]

SCHEMA_DIR = Path(workflow.current_basedir) / "../workflow/schemas"
GENOTYPES_TSV = Path(config["inputs"]["genotypes"])
GENOTYPES_DF = pl.read_csv(
    GENOTYPES_TSV,
    separator="\t",
    null_values="",
)
validate(
    GENOTYPES_DF,
    SCHEMA_DIR / "genotypes.schema.yaml",
    set_default=False,
)
check_genotypes_table(GENOTYPES_DF)
GENOTYPES = read_genotypes(GENOTYPES_TSV)
GENOTYPE_NAMES = [record["genotype"] for record in GENOTYPES]
REGION_FASTA_BY_GENOTYPE = {
    record["genotype"]: record["region_fasta"] for record in GENOTYPES
}

OUTDIR = Path(config["project"]["output_dir"])
SCRIPTS_DIR = Path(workflow.current_basedir) / "../scripts"
PROJECT_TITLE = config.get("project", {}).get("name", "Project")
VISUALIZATION_CONFIG = config.get("viewer", {})
SELECTED_MARKER_SETS = {}

SELECTED_MARKER_SETS_BY_SLUG = {
    slugify_marker_set_name(name): {
        "name": name,
        "tsv": path,
        "title": f"{PROJECT_TITLE} - {name}" if PROJECT_TITLE else name,
    }
    for name, path in SELECTED_MARKER_SETS.items()
}

if len(SELECTED_MARKER_SETS_BY_SLUG) != len(SELECTED_MARKER_SETS):
    raise ValueError("Several selected marker set names produce the same filename slug.")

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
DOTPLOT_PAIRS = resolve_dotplot_pairs(GENOTYPE_NAMES, config)
DOTPLOT_PAIR_IDS = [build_pair_id(sample_a, sample_b) for sample_a, sample_b in DOTPLOT_PAIRS]
DOTPLOT_PAFS = expand(
    DOTPLOT_PAF_DIR / "{pair_id}.paf",
    pair_id=DOTPLOT_PAIR_IDS,
)
DOTPLOT_FORMATTED = expand(
    DOTPLOT_FORMATTED_DIR / "{pair_id}.tsv",
    pair_id=DOTPLOT_PAIR_IDS,
)
# DOTPLOT_PDFS = expand(
#     DOTPLOT_PDF_DIR / "{pair_id}.pdf",
#     pair_id=DOTPLOT_PAIR_IDS,
# )
DOTPLOT_SVGS = expand(
    DOTPLOT_SVG_DIR / "{pair_id}.svg",
    pair_id=DOTPLOT_PAIR_IDS,
)
DOTPLOT_ONLY_SVGS = expand(
    DOTPLOT_ONLY_SVG_DIR / "{pair_id}.dotplot_only.svg",
    pair_id=DOTPLOT_PAIR_IDS,
)
DOTPLOT_GALLERY_HTML = DOTPLOT_COMBINED_DIR / "dotplots_gallery.html"
SUMMARY_STATS_JSON = SUMMARY_STATS_DIR / "summary_stats.json"
SUMMARY_STATS_TXT = SUMMARY_STATS_DIR / "summary_stats.txt"
DOTPLOT_MANIFEST = REGION_TRACK_DIR / "dotplots_manifest.json"
REGION_TRACK_HTML = REGION_TRACK_DIR / "region_tracks.html"
SELECTED_SNP_LONGS = expand(
    REGION_TRACK_DIR / "snp_positions_long.{marker_set}.tsv",
    marker_set=SELECTED_MARKER_SETS_BY_SLUG.keys(),
)
REGION_TRACK_SELECTED_HTMLS = expand(
    REGION_TRACK_DIR / "region_tracks.{marker_set}.html",
    marker_set=SELECTED_MARKER_SETS_BY_SLUG.keys(),
)
MASHTREE_MATRIX = MASH_DISTANCES_DIR / "mashtree.matrix.tsv"
MASHTREE_TREE = MASH_DISTANCES_DIR / "mashtree.dnd"
MASKED_BLOCK_N_STATS_TSV = BLOCK_STATS_DIR / "masked_block_n_stats.tsv"
GFF_TRACKS = get_gff_tracks(config, GENOTYPE_NAMES)
GFF_TRACK_FILES = get_gff_track_files(GFF_TRACKS)
GFF_TRACKS_JSON = REGION_TRACK_DIR / "gff_tracks.json"

NB_GENOTYPES = len(GENOTYPES)
te_lib_value = config.get("snps", {}).get("repeat_masking", {}).get("library", "")
TE_LIB = Path(te_lib_value) if te_lib_value else None
USE_MASKING = TE_LIB is not None

SNP_FILTER_GROUP_A, SNP_FILTER_GROUP_B, USE_SNP_GROUP_FILTERING = resolve_snp_filter_groups(
    genotypes=GENOTYPES,
)
