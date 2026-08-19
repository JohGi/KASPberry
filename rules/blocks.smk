rule rename_fasta_header:
    input:
        lambda wildcards: REGION_FASTA_BY_GENOTYPE[wildcards.sample]
    output:
        CLEAN_FASTA_DIR / "{sample}.fasta"
    benchmark:
        BENCHMARK_DIR / "rename_fasta_header" / "{sample}.tsv"
    log:
        stdout=LOG_DIR / "rename_fasta_header" / "{sample}.stdout",
        stderr=LOG_DIR / "rename_fasta_header" / "{sample}.stderr"
    shell:
        r"""
        mkdir -p "{CLEAN_FASTA_DIR}" "$(dirname "{log.stdout}")"
        bash "{SCRIPTS_DIR}/rename_fasta_header.sh" \
            --fasta "{input}" \
            --name "{wildcards.sample}" \
            --output "{output}" \
            1> "{log.stdout}" \
            2> "{log.stderr}"
        """

rule build_all_genomes_multifasta:
    input:
        CLEAN_FASTAS
    output:
        ALL_GENOMES_FASTA
    benchmark:
        BENCHMARK_DIR / "build_all_genomes_multifasta.tsv"
    log:
        stderr=LOG_DIR / "build_all_genomes_multifasta" / "build_all_genomes_multifasta.stderr"
    shell:
        r"""
        mkdir -p "{MULTIFASTA_DIR}" "$(dirname "{log.stderr}")"
        cat {input} > "{output}" \
            2> "{log.stderr}"
        """

rule run_sibeliaz:
    input:
        CLEAN_FASTAS
    output:
        SIBELIAZ_GFF
    benchmark:
        BENCHMARK_DIR / "run_sibeliaz.tsv"
    log:
        stdout=LOG_DIR / "run_sibeliaz" / "run_sibeliaz.stdout",
        stderr=LOG_DIR / "run_sibeliaz" / "run_sibeliaz.stderr"
    threads: 1
    params:
        outdir=str(SIBELIAZ_DIR),
        extra_options=config["advanced"]["sibeliaz"]["extra_options"]
    shell:
        r"""
        mkdir -p "{params.outdir}" "$(dirname "{log.stdout}")"
        sibeliaz \
            -n \
            -t {threads} \
            {params.extra_options} \
            -o "{params.outdir}" \
            {input} \
            1> "{log.stdout}" \
            2> "{log.stderr}"
        """

rule filter_sibeliaz_blocks:
    input:
        SIBELIAZ_GFF
    output:
        FILTERED_GFF
    benchmark:
        BENCHMARK_DIR / "filter_sibeliaz_blocks.tsv"
    log:
        stdout=LOG_DIR / "filter_sibeliaz_blocks" / "filter_sibeliaz_blocks.stdout",
        stderr=LOG_DIR / "filter_sibeliaz_blocks" / "filter_sibeliaz_blocks.stderr"
    params:
        nb_samples=NB_GENOTYPES,
        min_len=config["snps"]["min_block_length"]
    shell:
        r"""
        mkdir -p "$(dirname "{output}")" "$(dirname "{log.stdout}")"
        bash "{SCRIPTS_DIR}/filter_sibeliaz_blocks.sh" \
            --gff "{input}" \
            --nb_samples {params.nb_samples} \
            --min_len {params.min_len} \
            --output "{output}" \
            1> "{log.stdout}" \
            2> "{log.stderr}"
        """

rule map_block_coordinates:
    input:
        gff=FILTERED_GFF,
        genotypes_tsv=GENOTYPES_TSV
    output:
        BLOCK_COORDINATES_TSV
    benchmark:
        BENCHMARK_DIR / "map_block_coordinates.tsv"
    log:
        stdout=LOG_DIR / "map_block_coordinates" / "map_block_coordinates.stdout",
        stderr=LOG_DIR / "map_block_coordinates" / "map_block_coordinates.stderr"
    shell:
        r"""
        mkdir -p "$(dirname "{output}")" "$(dirname "{log.stdout}")"
        python3 "{SCRIPTS_DIR}/map_block_coordinates.py" \
            --gff "{input.gff}" \
            --genotypes-tsv "{input.genotypes_tsv}" \
            --output "{output}" \
            > "{log.stdout}" \
            2> "{log.stderr}"
        """

rule collect_blocks:
    input:
        FILTERED_GFF
    output:
        BLOCK_LIST
    benchmark:
        BENCHMARK_DIR / "collect_blocks.tsv"
    log:
        stderr=LOG_DIR / "collect_blocks" / "collect_blocks.stderr",
        stdout=LOG_DIR / "collect_blocks" / "collect_blocks.stdout"
    shell:
        r"""
        mkdir -p "{FILTERED_BLOCKS_DIR}" "$(dirname "{log.stdout}")"
        bash "{SCRIPTS_DIR}/extract_block_ids.sh" \
            "{input}" \
            "{output}" \
            1> "{log.stdout}" \
            2> "{log.stderr}"
        """

rule extract_all_blocks_fasta:
    input:
        fasta=ALL_GENOMES_FASTA,
        gff=FILTERED_GFF
    output:
        ALL_BLOCKS_RAW_FASTA
    benchmark:
        BENCHMARK_DIR / "extract_all_blocks_fasta.tsv"
    log:
        stdout=LOG_DIR / "extract_all_blocks_fasta" / "extract_all_blocks_fasta.stdout",
        stderr=LOG_DIR / "extract_all_blocks_fasta" / "extract_all_blocks_fasta.stderr"
    shell:
        r"""
        mkdir -p "{BLOCK_FASTA_DIR}" "$(dirname "{log.stdout}")"
        bash "{SCRIPTS_DIR}/extract_all_blocks_fasta.sh" \
            --fasta "{input.fasta}" \
            --gff "{input.gff}" \
            --output "{output}" \
            1> "{log.stdout}" \
            2> "{log.stderr}"
        """


checkpoint split_block_fastas:
    input:
        raw_fasta=ALL_BLOCKS_RAW_FASTA,
        block_list=BLOCK_LIST
    output:
        directory(BLOCK_FASTA_SPLIT_DIR)
    benchmark:
        BENCHMARK_DIR / "split_block_fastas.tsv"
    log:
        stdout=LOG_DIR / "split_block_fastas" / "split_block_fastas.stdout",
        stderr=LOG_DIR / "split_block_fastas" / "split_block_fastas.stderr"
    shell:
        r"""
        mkdir -p "{BLOCK_FASTA_SPLIT_DIR}" "$(dirname "{log.stdout}")"
        python3 "{SCRIPTS_DIR}/split_block_fastas.py" \
            --input "{input.raw_fasta}" \
            --block-list "{input.block_list}" \
            --outdir "{output}" \
            1> "{log.stdout}" \
            2> "{log.stderr}"
        """

checkpoint split_block_list_into_chunks:
    input:
        BLOCK_LIST
    output:
        chunk_dir=directory(MASK_CHUNK_DIR)
    benchmark:
        BENCHMARK_DIR / "split_block_list_into_chunks.tsv"
    log:
        stdout=LOG_DIR / "split_block_list_into_chunks" / "split_block_list_into_chunks.stdout",
        stderr=LOG_DIR / "split_block_list_into_chunks" / "split_block_list_into_chunks.stderr"
    params:
        chunk_size=config["advanced"]["batching"]["blocks_per_chunk"]
    shell:
        r"""
        mkdir -p "{FILTERED_BLOCKS_DIR}" "$(dirname "{log.stdout}")"
        python3 "{SCRIPTS_DIR}/split_block_list_into_chunks.py" \
            --input "{input}" \
            --output-dir "{output.chunk_dir}" \
            --chunk-size {params.chunk_size} \
            1> "{log.stdout}" \
            2> "{log.stderr}"
        """
