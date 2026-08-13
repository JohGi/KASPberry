POLYMARKER_INPUT_DIR = OUTDIR / "15_polymarker_inputs"

POLYMARKER_GENOME_FASTA = (
    POLYMARKER_INPUT_DIR / "{genotype}" / "genome.fasta"
)

POLYMARKER_ALIASES_TSV = (
    POLYMARKER_INPUT_DIR / "{genotype}" / "chromosome_aliases.tsv"
)


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


def get_kasp_genome_fasta(wildcards):
    """Return the whole-genome FASTA for one KASP QC genotype."""
    return KASP_GENOME_FASTA_BY_GENOTYPE[wildcards.genotype]


def get_chromosomes_input(_wildcards=None):
    """Return chromosomes.tsv when provided."""
    chromosomes_path = config.get("inputs", {}).get("chromosomes")

    if chromosomes_path:
        return chromosomes_path

    return []


rule prepare_polymarker_genome:
    input:
        genotypes=GENOTYPES_TSV,
        genome=get_kasp_genome_fasta,
        chromosomes=get_chromosomes_input,
    output:
        genome=POLYMARKER_GENOME_FASTA,
        aliases=POLYMARKER_ALIASES_TSV,
    params:
        polymarker_genomes=config["kasp"]["polymarker_genomes"],
    log:
        stdout=LOG_DIR / "prepare_polymarker_genome" / "{genotype}.stdout",
        stderr=LOG_DIR / "prepare_polymarker_genome" / "{genotype}.stderr",
    benchmark:
        BENCHMARK_DIR / "prepare_polymarker_genome" / "{genotype}.tsv"
    script:
        "python3 "{SCRIPTS_DIR}/prepare_polymarker_genome.py ..."
