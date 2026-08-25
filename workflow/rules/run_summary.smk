import os


# The KASPberry CLI sets this explicitly; Direct Snakefile invocations default to "snps"
RUN_SUMMARY_MODE = os.environ.get("KASPBERRY_MODE", "snps")


def get_run_summary_mode():
    """Return the requested workflow mode."""
    if RUN_SUMMARY_MODE not in {"snps", "kasp"}:
        raise ValueError(
            "Invalid KASPBERRY_MODE. Expected 'snps' or 'kasp', "
            f"got {RUN_SUMMARY_MODE!r}."
        )

    return RUN_SUMMARY_MODE


def get_run_summary_snp_summary(_wildcards):
    """Return the aggregated SNP summary for the current mode."""
    mode = get_run_summary_mode()

    if mode == "snps":
        return SNP_SUMMARY_TSV

    return KASP_SUMMARY_TSV


def get_run_summary_assay_summary(_wildcards):
    """Return the assay summary only for KASP runs."""
    if get_run_summary_mode() == "kasp":
        return ASSAY_SUMMARY_TSV

    return []


rule write_run_summary:
    conda:
        "../envs/workflow-runtime.yaml"
    input:
        genotypes=GENOTYPES_TSV,
        block_coords=BLOCK_COORDINATES_TSV,
        snp_summary=get_run_summary_snp_summary,
        assay_summary=get_run_summary_assay_summary,
        clean_fastas=CLEAN_FASTAS,
        masked_block_n_stats=MASKED_BLOCK_N_STATS_TSV,
    output:
        json=RUN_SUMMARY_JSON,
        txt=RUN_SUMMARY_TXT
    params:
        mode=lambda wildcards: get_run_summary_mode(),
        repeat_library=str(TE_LIB) if TE_LIB is not None else "",
    benchmark:
        BENCHMARK_DIR / "write_run_summary.tsv"
    log:
        stdout=LOG_DIR / "write_run_summary" / "write_run_summary.stdout",
        stderr=LOG_DIR / "write_run_summary" / "write_run_summary.stderr"
    shell:
        r"""
        mkdir -p \
            "$(dirname "{output.json}")" \
            "$(dirname "{output.txt}")" \
            "$(dirname "{log.stdout}")"

        assay_args=()

        if [[ "{params.mode}" == "kasp" ]]; then
            assay_args=(
                --assay-summary "{input.assay_summary}"
            )
        fi

        repeat_masking_args=()

        if [[ -n "{params.repeat_library}" ]]; then
            repeat_masking_args=(
                --repeat-masking-library "{params.repeat_library}"
            )
        fi

        python3 "{SCRIPTS_DIR}/write_run_summary.py" \
            --mode "{params.mode}" \
            --genotypes "{input.genotypes}" \
            --block-coords "{input.block_coords}" \
            --snp-summary "{input.snp_summary}" \
            "${{assay_args[@]}}" \
            --repeat-masking \
            "${{repeat_masking_args[@]}}" \
            --clean-fastas {input.clean_fastas:q} \
            --masked-block-n-stats "{input.masked_block_n_stats}" \
            --json-output "{output.json}" \
            --txt-output "{output.txt}" \
            1> "{log.stdout}" \
            2> "{log.stderr}"
        """
