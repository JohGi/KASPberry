rule compute_masked_block_n_stats:
    input:
        masked_align_done=get_align_chunk_sentinels,
        unmasked_align_done=get_unmasked_align_chunk_sentinels
    output:
        MASKED_BLOCK_N_STATS_TSV
    benchmark:
        BENCHMARK_DIR / "compute_masked_block_n_stats.tsv"
    log:
        stdout=LOG_DIR / "compute_masked_block_n_stats" / "compute_masked_block_n_stats.stdout",
        stderr=LOG_DIR / "compute_masked_block_n_stats" / "compute_masked_block_n_stats.stderr"
    shell:
        r"""
        mkdir -p "{BLOCK_STATS_DIR}" "$(dirname "{log.stdout}")"
        bash "{SCRIPTS_DIR}/write_masked_block_n_stats.sh" \
            --masked-dir "{ALIGN_DIR}" \
            --unmasked-dir "{UNMASKED_ALIGN_DIR}" \
            --output "{output}" \
            > "{log.stdout}" \
            2> "{log.stderr}"
        """

rule align_unmasked_block_chunk:
    input:
        chunk_list=MASK_CHUNK_DIR / "{chunk_id}.list",
        split_block_dir=get_split_block_dir
    output:
        UNMASKED_ALIGN_DIR / "{chunk_id}.done"
    benchmark:
        BENCHMARK_DIR / "align_unmasked_block_chunk" / "{chunk_id}.tsv"
    log:
        stdout=LOG_DIR / "align_unmasked_block_chunk" / "{chunk_id}.stdout",
        stderr=LOG_DIR / "align_unmasked_block_chunk" / "{chunk_id}.stderr"
    threads: 1
    params:
        fasta_dir=str(BLOCK_FASTA_SPLIT_DIR),
        outdir=str(UNMASKED_ALIGN_DIR),
        fasta_suffix=".fasta",
        extra_options=" ".join(config.get("advanced", {}).get("alignment", {}).get("mafft_options", []))
    shell:
        r"""
        mkdir -p "{params.outdir}" "$(dirname "{log.stdout}")"
        bash "{SCRIPTS_DIR}/align_block_chunk.sh" \
            --chunk-list "{input.chunk_list}" \
            --fasta-dir "{params.fasta_dir}" \
            --outdir "{params.outdir}" \
            --threads {threads} \
            --fasta-suffix "{params.fasta_suffix}" \
            --mafft-extra-options "{params.extra_options}" \
            1> "{log.stdout}" \
            2> "{log.stderr}"
        touch "{output}"
        """

rule compute_kimura2p_distmat_chunk:
    input:
        chunk_list=MASK_CHUNK_DIR / "{chunk_id}.list",
        unmasked_align_done=UNMASKED_ALIGN_DIR / "{chunk_id}.done"
    output:
        KIMURA2P_DISTMAT_CHUNK_DIR / "{chunk_id}.done"
    benchmark:
        BENCHMARK_DIR / "compute_kimura2p_distmat_chunk" / "{chunk_id}.tsv"
    log:
        stdout=LOG_DIR / "compute_kimura2p_distmat_chunk" / "{chunk_id}.stdout",
        stderr=LOG_DIR / "compute_kimura2p_distmat_chunk" / "{chunk_id}.stderr"
    params:
        aln_dir=str(UNMASKED_ALIGN_DIR),
        outdir=str(KIMURA2P_DISTMAT_MATRIX_DIR)
    shell:
        r"""
        mkdir -p "{params.outdir}" "{KIMURA2P_DISTMAT_CHUNK_DIR}" "$(dirname "{log.stdout}")" && \
        bash "{SCRIPTS_DIR}/compute_kimura2p_distmat_chunk.sh" \
            --chunk-list "{input.chunk_list}" \
            --aln-dir "{params.aln_dir}" \
            --outdir "{params.outdir}" \
            1> "{log.stdout}" \
            2> "{log.stderr}" && \
        touch "{output}"
        """
