rule write_run_summary:
    input:
        block_coords=BLOCK_COORDINATES_TSV,
        snp_positions=SNP_POS_WIDE_TSV,
        clean_fastas=CLEAN_FASTAS
    output:
        json=RUN_SUMMARY_JSON,
        txt=RUN_SUMMARY_TXT
    benchmark:
        BENCHMARK_DIR / "write_run_summary.tsv"
    log:
        stdout=LOG_DIR / "write_run_summary" / "write_run_summary.stdout",
        stderr=LOG_DIR / "write_run_summary" / "write_run_summary.stderr"
    shell:
        r"""
        mkdir -p "$(dirname "{output.txt}")" "$(dirname "{log.stdout}")"
        python3 "{SCRIPTS_DIR}/write_run_summary.py" \
            --block-coords "{input.block_coords}" \
            --snp-positions "{input.snp_positions}" \
            --clean-fastas {input.clean_fastas} \
            --json-output "{output.json}" \
            --txt-output "{output.txt}" \
            > "{log.stdout}" \
            2> "{log.stderr}"
        """
