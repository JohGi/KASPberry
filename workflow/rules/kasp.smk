def get_kasp_genome_fasta(wildcards):
    """Return the whole-genome FASTA for one KASP QC genotype."""
    return KASP_GENOME_FASTA_BY_GENOTYPE[wildcards.genotype]


def get_chromosomes_input(_wildcards=None):
    """Return chromosomes.tsv when provided."""
    chromosomes_path = config["inputs"]["chromosomes"]

    if chromosomes_path:
        return chromosomes_path

    return []


def get_polymarker_fai(wildcards):
    return f"{POLYMARKER_INPUT_DIR / wildcards.genotype / 'genome.fasta'}.fai"


rule prepare_polymarker_genome:
    conda:
        "../envs/kasp-preparation.yaml"
    input:
        genome=get_kasp_genome_fasta,
        chromosomes=get_chromosomes_input,
    output:
        genome=POLYMARKER_INPUT_DIR / "{genotype}/genome.fasta",
        aliases=POLYMARKER_INPUT_DIR / "{genotype}/chromosome_aliases.tsv",
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
    conda:
        "../envs/kasp-preparation.yaml"
    input:
        genome=POLYMARKER_INPUT_DIR / "{genotype}/genome.fasta",
    output:
        fai=POLYMARKER_INPUT_DIR / "{genotype}/genome.fasta.fai",
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
    conda:
        "../envs/kasp-preparation.yaml"
    input:
        genome=POLYMARKER_INPUT_DIR / "{genotype}/genome.fasta",
    output:
        done=POLYMARKER_INPUT_DIR / "{genotype}/genome.fasta.blastdb.done",
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
    conda:
        "../envs/kasp-preparation.yaml"
    input:
        snps=SNP_POS_LONG_TSV,
        snp_status=DIAGNOSTIC_STATUS_TSV,
        reference=CLEAN_FASTA_DIR / "{genotype}.fasta",
        aliases=POLYMARKER_INPUT_DIR / "{genotype}/chromosome_aliases.tsv",
    output:
        marker_list=POLYMARKER_INPUT_DIR / "{genotype}/marker_list.csv",
    params:
        flank=POLYMARKER_CONTEXT_FLANK,
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
            --snp-status "{input.snp_status}" \
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
        genomes=lambda wildcards: config["kasp"]["polymarker_genomes"],
    container:
            "docker://ghcr.io/johgi/kaspberry-polymarker:1.3.3-k1"
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
            --arm_selection first_two \
            --output "{params.outdir}" \
            --genomes_count "{params.genomes}" \
            --primers_to_order \
            1> "{log.stdout}" \
            2> "{log.stderr}"
        """


rule process_polymarker_outputs:
    conda:
        "../envs/kasp-preparation.yaml"
    input:
        marker_lists=expand(
            POLYMARKER_INPUT_DIR / "{genotype}/marker_list.csv",
            genotype=KASP_GENOTYPE_NAMES,
        ),
        primers_to_order=expand(
            POLYMARKER_DIR / "{genotype}/primers_to_order.csv",
            genotype=KASP_GENOTYPE_NAMES,
        ),
    output:
        snp_status=POLYMARKER_DESIGN_STATUS_TSV,
        snp_status_by_genotype=POLYMARKER_DESIGN_STATUS_BY_GENOTYPE_TSV,
        assays=POLYMARKER_ASSAYS_TSV,
        mfe_canonical_pairs=POLYMARKER_SUMMARY_DIR / "primers_canonical_pairs.tsv",
        mfe_noncanonical_pairs=POLYMARKER_SUMMARY_DIR / "primers_noncanonical_pairs.tsv",
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
            --snp-status "{output.snp_status}" \
            --snp-status-by-genotype "{output.snp_status_by_genotype}" \
            --assays "{output.assays}" \
            --mfe-canonical-pairs "{output.mfe_canonical_pairs}" \
            --mfe-noncanonical-pairs "{output.mfe_noncanonical_pairs}" \
            --mfe-with-tails "{output.mfe_with_tails}" \
            1> "{log.stdout}" \
            2> "{log.stderr}"
        """
