def resolve_snp_filter_groups(
    genotypes: list[dict[str, str]],
) -> tuple[list[str], list[str]]:
    """Resolve SNP filtering groups from the genotype table."""
    groups = sorted({row["group"] for row in genotypes if row.get("group")})

    if not groups:
        return [], []

    group_a = [row["genotype"] for row in genotypes if row.get("group") == groups[0]]
    group_b = [row["genotype"] for row in genotypes if row.get("group") == groups[1]]

    return group_a, group_b


SNP_FILTER_GROUP_A, SNP_FILTER_GROUP_B = resolve_snp_filter_groups(
    genotypes=GENOTYPES,
)


rule mask_block_chunk:
    input:
        fasta_dir=get_split_block_dir,
        chunk_list=MASK_CHUNK_DIR / "{chunk_id}.list"
    output:
        MASKED_CHUNK_DIR / "{chunk_id}.done"
    benchmark:
        BENCHMARK_DIR / "mask_block_chunk" / "{chunk_id}.tsv"
    log:
        stdout=LOG_DIR / "mask_block_chunk" / "{chunk_id}.stdout",
        stderr=LOG_DIR / "mask_block_chunk" / "{chunk_id}.stderr"
    threads: 1
    params:
        te_lib=TE_LIB,
        outdir=str(MASKED_DIR)
    shell:
        r"""
        mkdir -p "{MASKED_CHUNK_DIR}" "$(dirname "{log.stdout}")"
        bash "{SCRIPTS_DIR}/mask_block_chunk.sh" \
            --chunk-list "{input.chunk_list}" \
            --fasta-dir "{input.fasta_dir}" \
            --te-lib "{params.te_lib}" \
            --outdir "{params.outdir}" \
            --threads {threads} \
            1> "{log.stdout}" \
            2> "{log.stderr}"
        touch "{output}"
        """


rule align_block_chunk:
    input:
        get_alignment_inputs
    output:
        ALIGN_DIR / "{chunk_id}.done"
    benchmark:
        BENCHMARK_DIR / "align_block_chunk" / "{chunk_id}.tsv"
    log:
        stdout=LOG_DIR / "align_block_chunk" / "{chunk_id}.stdout",
        stderr=LOG_DIR / "align_block_chunk" / "{chunk_id}.stderr"
    threads: 1
    params:
        fasta_dir=str(ALIGNMENT_FASTA_DIR),
        outdir=str(ALIGN_DIR),
        extra_options=" ".join(config["advanced"]["alignment"]["mafft_options"])
    shell:
        r"""
        mkdir -p "{ALIGN_DIR}" "$(dirname "{log.stdout}")"
        bash "{SCRIPTS_DIR}/align_block_chunk.sh" \
            --chunk-list "{input[0]}" \
            --fasta-dir "{params.fasta_dir}" \
            --outdir "{params.outdir}" \
            --threads {threads} \
            --fasta-suffix .fasta.masked \
            --mafft-extra-options "{params.extra_options}" \
            1> "{log.stdout}" \
            2> "{log.stderr}" && \
        touch "{output}"
        """

rule detect_snps:
    input:
        align_sentinels=get_align_chunk_sentinels
    output:
        SNP_VCF
    benchmark:
        BENCHMARK_DIR / "detect_snps.tsv"
    log:
        stdout=LOG_DIR / "detect_snps" / "detect_snps.stdout",
        stderr=LOG_DIR / "detect_snps" / "detect_snps.stderr"
    params:
        min_flank=config["snps"]["min_snp_flank"]
    shell:
        r"""
        mkdir -p "{SNP_DIR}" "$(dirname "{log.stdout}")"
        bash "{SCRIPTS_DIR}/detect_snps.sh" \
            --align-dir "{ALIGN_DIR}" \
            --output "{output}" \
            --min-flank "{params.min_flank}" \
            1> "{log.stdout}" \
            2> "{log.stderr}"
        """

rule write_snp_group_files:
    output:
        group_a=GROUP_A_LIST,
        group_b=GROUP_B_LIST
    benchmark:
        BENCHMARK_DIR / "write_snp_group_files.tsv"
    run:
        Path(output.group_a).parent.mkdir(parents=True, exist_ok=True)
        Path(output.group_a).write_text(
            "\n".join(SNP_FILTER_GROUP_A) + ("\n" if SNP_FILTER_GROUP_A else ""),
            encoding="utf-8",
        )
        Path(output.group_b).write_text(
            "\n".join(SNP_FILTER_GROUP_B) + ("\n" if SNP_FILTER_GROUP_B else ""),
            encoding="utf-8",
        )

rule filter_snps_by_groups:
    input:
        snps=SNP_VCF,
        group_a=GROUP_A_LIST,
        group_b=GROUP_B_LIST
    output:
        status=DIAGNOSTIC_STATUS_TSV
    benchmark:
        BENCHMARK_DIR / "filter_snps_by_groups.tsv"
    log:
        stdout=LOG_DIR / "filter_snps_by_groups" / "filter_snps_by_groups.stdout",
        stderr=LOG_DIR / "filter_snps_by_groups" / "filter_snps_by_groups.stderr"
    shell:
        r"""
        mkdir -p "{SNP_WORK_DIR}" "$(dirname "{log.stdout}")"
        python3 "{SCRIPTS_DIR}/filter_snps_by_groups.py" \
            --input "{input.snps}" \
            --output "{output.status}" \
            --group-a-file "{input.group_a}" \
            --group-b-file "{input.group_b}" \
            1> "{log.stdout}" \
            2> "{log.stderr}"
        """
