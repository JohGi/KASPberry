rule map_snp_coordinates:
    input:
        vcf=get_final_snp_output(),
        block_coords=BLOCK_COORDINATES_TSV,
        samples_tsv=GENOTYPES_TSV
    output:
        long=SNP_POS_LONG_TSV,
        wide=SNP_POS_WIDE_TSV
    benchmark:
        BENCHMARK_DIR / "map_snp_coordinates.tsv"
    log:
        stdout=LOG_DIR / "map_snp_coordinates" / "map_snp_coordinates.stdout",
        stderr=LOG_DIR / "map_snp_coordinates" / "map_snp_coordinates.stderr"
    shell:
        r"""
        mkdir -p "$(dirname {output.long})" "$(dirname "{log.stdout}")"
        python3 "{SCRIPTS_DIR}/map_snp_coordinates.py" \
            --vcf "{input.vcf}" \
            --block-coords "{input.block_coords}" \
            --samples-tsv "{input.samples_tsv}" \
            --align-dir "{ALIGN_DIR}" \
            --long-output "{output.long}" \
            --wide-output "{output.wide}" \
            > "{log.stdout}" \
            2> "{log.stderr}"
        """
