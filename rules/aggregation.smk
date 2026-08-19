rule aggregate_snp_results:
    input:
        vcf=SNP_VCF,
        positions=SNP_POS_WIDE_TSV,
        diagnostic_status=DIAGNOSTIC_STATUS_TSV,
    output:
        summary=SNP_SUMMARY_TSV,
        retained_vcf=DIAGNOSTIC_SNPS_VCF,
    shell:
        """
        mkdir -p "$(dirname \"{output.summary}\")" "$(dirname \"{output.retained_vcf}\")"
        python3 "{SCRIPTS_DIR}/aggregate_results.py" \
            --mode snps \
            --snps-vcf "{input.vcf}" \
            --snp-positions "{input.positions}" \
            --diagnostic-status "{input.diagnostic_status}" \
            --output "{output.summary}" \
            --retained-vcf "{output.retained_vcf}"
        """


rule aggregate_kasp_results:
    input:
        vcf=SNP_VCF,
        positions=SNP_POS_WIDE_TSV,
        diagnostic_status=DIAGNOSTIC_STATUS_TSV,
        design_status=POLYMARKER_DESIGN_STATUS_TSV,
        validation_status=IN_SILICO_VALIDATION_STATUS_TSV,
        assays=POLYMARKER_ASSAYS_TSV,
        assay_status=IN_SILICO_ASSAY_STATUS_TSV,
    output:
        summary=KASP_SUMMARY_TSV,
        retained_vcf=CANDIDATE_SNPS_VCF,
        assay_summary=ASSAY_SUMMARY_TSV,
        candidate_assays=CANDIDATE_ASSAYS_TSV,
        primers_to_order=PRIMERS_TO_ORDER_TSV,
    shell:
        """
        mkdir -p "$(dirname \"{output.summary}\")" "$(dirname \"{output.retained_vcf}\")" "$(dirname \"{output.assay_summary}\")"
        python3 "{SCRIPTS_DIR}/aggregate_results.py" \
            --mode kasp \
            --snps-vcf "{input.vcf}" \
            --snp-positions "{input.positions}" \
            --diagnostic-status "{input.diagnostic_status}" \
            --design-status "{input.design_status}" \
            --validation-status "{input.validation_status}" \
            --assays "{input.assays}" \
            --assay-status "{input.assay_status}" \
            --output "{output.summary}" \
            --retained-vcf "{output.retained_vcf}" \
            --assay-summary "{output.assay_summary}" \
            --candidate-assays "{output.candidate_assays}" \
            --primers-to-order "{output.primers_to_order}"
        """
