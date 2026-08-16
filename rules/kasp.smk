POLYMARKER_INPUT_DIR = OUTDIR / "15_polymarker_inputs"
POLYMARKER_DIR = OUTDIR / "16_polymarker"
POLYMARKER_SUMMARY_DIR = OUTDIR / "17_polymarker_summary"

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


# POLYMARKER_MARKER_LIST = POLYMARKER_INPUT_DIR / "marker_list.csv"
# POLYMARKER_BLASTDB_DONE = POLYMARKER_INPUT_DIR / "{genotype}" / "blastdb.done"

# POLYMARKER_PRIMERS_CSV = POLYMARKER_DIR / "{genotype}" / "primers.csv"
# POLYMARKER_PRIMERS_TO_ORDER = POLYMARKER_DIR / "{genotype}" / "primers_to_order.csv"
# POLYMARKER_STATUS_TXT = POLYMARKER_DIR / "{genotype}" / "status.txt"

# POLYMARKER_SNP_STATUS_TSV = POLYMARKER_SUMMARY_DIR / "polymarker_snp_status.tsv"
# POLYMARKER_SNP_STATUS_LONG_TSV = POLYMARKER_SUMMARY_DIR / "polymarker_snp_status_long.tsv"
# POLYMARKER_ASSAYS_TSV = POLYMARKER_SUMMARY_DIR / "polymarker_assays.tsv"
# MFE_PRIMER_PAIRS_TSV = POLYMARKER_SUMMARY_DIR / "primers.tsv"
# MFE_CONTROL_PAIRS_TSV = POLYMARKER_SUMMARY_DIR / "primers_control_pairs.tsv"
# MFE_PRIMERS_WITH_TAILS_FASTA = POLYMARKER_SUMMARY_DIR / "primers_with_tails.fasta"


def get_kasp_genome_fasta(wildcards):
    """Return the whole-genome FASTA for one KASP QC genotype."""
    return KASP_GENOME_FASTA_BY_GENOTYPE[wildcards.genotype]


def get_chromosomes_input(_wildcards=None):
    """Return chromosomes.tsv when provided."""
    chromosomes_path = config.get("inputs", {}).get("chromosomes")

    if chromosomes_path:
        return chromosomes_path

    return []


def get_polymarker_fai(wildcards):
    return f"{POLYMARKER_INPUT_DIR / wildcards.genotype / 'genome.fasta'}.fai"


rule prepare_polymarker_genome:
    input:
        genome=get_kasp_genome_fasta,
        chromosomes=get_chromosomes_input,
    output:
        genome=OUTDIR / "15_polymarker_inputs/{genotype}/genome.fasta",
        aliases=OUTDIR / "15_polymarker_inputs/{genotype}/chromosome_aliases.tsv",
    params:
        source_seq=lambda wildcards: next(
            record["source_seq"]
            for record in KASP_GENOTYPES
            if record["genotype"] == wildcards.genotype
        ),
    log:
        stdout=LOG_DIR / "prepare_polymarker_genome/{genotype}.stdout",
        stderr=LOG_DIR / "prepare_polymarker_genome/{genotype}.stderr",
    benchmark:
        BENCHMARK_DIR / "prepare_polymarker_genome/{genotype}.tsv"
    shell:
        r"""
        mkdir -p "$(dirname "{output.genome}")" "$(dirname "{log.stdout}")"

        python3 "{SCRIPTS_DIR}/prepare_polymarker_genome.py" \
            --genotype "{wildcards.genotype}" \
            --source-seq "{params.source_seq}" \
            --genome "{input.genome}" \
            --chromosomes "{input.chromosomes}" \
            --out-fasta "{output.genome}" \
            --out-aliases "{output.aliases}" \
            > "{log.stdout}" \
            2> "{log.stderr}"
        """


rule index_polymarker_genome:
    input:
        genome=OUTDIR / "15_polymarker_inputs/{genotype}/genome.fasta",
    output:
        fai=OUTDIR / "15_polymarker_inputs/{genotype}/genome.fasta.fai",
    log:
        stdout=LOG_DIR / "index_polymarker_genome/{genotype}.stdout",
        stderr=LOG_DIR / "index_polymarker_genome/{genotype}.stderr",
    benchmark:
        BENCHMARK_DIR / "index_polymarker_genome/{genotype}.tsv"
    shell:
        r"""
        mkdir -p "$(dirname "{log.stdout}")"

        samtools faidx "{input.genome}" \
            1> "{log.stdout}" \
            2> "{log.stderr}"
        """


rule build_polymarker_blastdb:
    input:
        genome=OUTDIR / "15_polymarker_inputs/{genotype}/genome.fasta",
    output:
        done=OUTDIR / "15_polymarker_inputs/{genotype}/genome.fasta.blastdb.done",
    log:
        stdout=LOG_DIR / "build_polymarker_blastdb/{genotype}.stdout",
        stderr=LOG_DIR / "build_polymarker_blastdb/{genotype}.stderr",
    benchmark:
        BENCHMARK_DIR / "build_polymarker_blastdb/{genotype}.tsv"
    shell:
        r"""
        mkdir -p "$(dirname "{log.stdout}")"

        makeblastdb \
            -in "{input.genome}" \
            -dbtype nucl \
            -parse_seqids \
            1> "{log.stdout}" \
            2> "{log.stderr}"

        touch "{output.done}"
        """


rule prepare_polymarker_markers:
    input:
        snps=SNP_POS_LONG_TSV,
        reference=CLEAN_FASTA_DIR / "{genotype}.fasta",
        aliases=OUTDIR / "15_polymarker_inputs/{genotype}/chromosome_aliases.tsv",
    output:
        marker_list=OUTDIR / "15_polymarker_inputs/{genotype}/marker_list.csv",
    params:
        flank=config["snps"]["min_snp_flank"],
    log:
        stdout=LOG_DIR / "prepare_polymarker_markers/{genotype}.stdout",
        stderr=LOG_DIR / "prepare_polymarker_markers/{genotype}.stderr",
    benchmark:
        BENCHMARK_DIR / "prepare_polymarker_markers/{genotype}.tsv"
    shell:
        r"""
        mkdir -p "$(dirname "{output.marker_list}")" "$(dirname "{log.stdout}")"

        python3 "{SCRIPTS_DIR}/prepare_polymarker_markers.py" \
            --snps "{input.snps}" \
            --reference "{input.reference}" \
            --aliases "{input.aliases}" \
            --flank {params.flank} \
            --out "{output.marker_list}" \
            > "{log.stdout}" \
            2> "{log.stderr}"
        """


rule run_polymarker:
    input:
        genome=POLYMARKER_INPUT_DIR / "{genotype}/genome.fasta",
        fai=POLYMARKER_INPUT_DIR / "{genotype}/genome.fasta.fai",
        blastdb_done=POLYMARKER_INPUT_DIR / "{genotype}/genome.fasta.blastdb.done",
        marker_list=POLYMARKER_INPUT_DIR / "{genotype}/marker_list.csv",
    output:
        primers=POLYMARKER_DIR / "{genotype}/primers.csv",
        primers_to_order=POLYMARKER_DIR / "{genotype}/primers_to_order.csv",
    params:
        outdir=lambda wildcards: POLYMARKER_DIR / wildcards.genotype,
        genomes=config["kasp"]["polymarker_genomes"],
    container:
        "containers/bio-polyploid-tools.sif"
    threads: 1
    log:
        stdout=LOG_DIR / "run_polymarker/{genotype}.stdout",
        stderr=LOG_DIR / "run_polymarker/{genotype}.stderr",
    benchmark:
        BENCHMARK_DIR / "run_polymarker/{genotype}.tsv"
    shell:
        r"""
        mkdir -p "{params.outdir}" "$(dirname "{log.stdout}")"

        polymarker.rb \
            --contigs "{input.genome}" \
            --marker_list "{input.marker_list}" \
            --arm_selection arm_selection_first_two \
            --output "{params.outdir}" \
            --genomes "{params.genomes}" \
            --primers_to_order \
            1> "{log.stdout}" \
            2> "{log.stderr}"
        """


rule process_polymarker_outputs:
    input:
        marker_lists=expand(
            POLYMARKER_INPUT_DIR / "{genotype}/marker_list.csv",
            genotype=KASP_GENOTYPE_NAMES,
        ),
        primers_to_order=expand(
            POLYMARKER_DIR / "{genotype}/primers_to_order.csv",
            genotype=KASP_GENOTYPE_NAMES,
        ),
        primers_csv=expand(
            POLYMARKER_DIR / "{genotype}/primers.csv",
            genotype=KASP_GENOTYPE_NAMES,
        ),
    output:
        snp_status=POLYMARKER_SUMMARY_DIR / "polymarker_snp_status.tsv",
        snp_status_long=POLYMARKER_SUMMARY_DIR / "polymarker_snp_status_long.tsv",
        assays=POLYMARKER_SUMMARY_DIR / "polymarker_assays.tsv",
        mfe_pairs=POLYMARKER_SUMMARY_DIR / "primers.tsv",
        mfe_control_pairs=POLYMARKER_SUMMARY_DIR / "primers_control_pairs.tsv",
        mfe_with_tails=POLYMARKER_SUMMARY_DIR / "primers_with_tails.fasta",
    log:
        stdout=LOG_DIR / "process_polymarker_outputs.stdout",
        stderr=LOG_DIR / "process_polymarker_outputs.stderr",
    benchmark:
        BENCHMARK_DIR / "process_polymarker_outputs.tsv"
    shell:
        r"""
        mkdir -p "{POLYMARKER_SUMMARY_DIR}" "$(dirname "{log.stdout}")"

        python3 "{SCRIPTS_DIR}/process_polymarker_outputs.py" \
            --marker-lists {input.marker_lists:q} \
            --primers-to-order {input.primers_to_order:q} \
            --primers-csv {input.primers_csv:q} \
            --snp-status "{output.snp_status}" \
            --snp-status-long "{output.snp_status_long}" \
            --assays "{output.assays}" \
            --mfe-pairs "{output.mfe_pairs}" \
            --mfe-control-pairs "{output.mfe_control_pairs}" \
            --mfe-with-tails "{output.mfe_with_tails}" \
            1> "{log.stdout}" \
            2> "{log.stderr}"
        """
