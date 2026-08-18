rule index_mfeprimer_genome:
    input:
        genome=POLYMARKER_INPUT_DIR / "{genotype}/genome.fasta",
    output:
        db=IN_SILICO_DIR / "indexes/{genotype}/genome.fasta",
        done=IN_SILICO_DIR / "indexes/{genotype}/index.done",
    log:
        stdout=LOG_DIR / "index_mfeprimer_genome/{genotype}.stdout",
        stderr=LOG_DIR / "index_mfeprimer_genome/{genotype}.stderr",
    benchmark:
        BENCHMARK_DIR / "index_mfeprimer_genome/{genotype}.tsv"
    shell:
        r"""
        mkdir -p "$(dirname "{output.db}")" "$(dirname "{log.stdout}")"

        ln -srf "{input.genome}" "{output.db}"

        "{workflow.basedir}/dependencies/bin/mfeprimer-4.2.3-linux-amd64" index \
            -i "{output.db}" \
            1> "{log.stdout}" \
            2> "{log.stderr}"

        touch "{output.done}"
        """


rule run_mfeprimer_specificity:
    wildcard_constraints:
        pair_set="canonical|noncanonical"

    input:
        primers=lambda wildcards: {
            "canonical":
                POLYMARKER_SUMMARY_DIR / "primers_canonical_pairs.tsv",
            "noncanonical":
                POLYMARKER_SUMMARY_DIR / "primers_noncanonical_pairs.tsv",
        }[wildcards.pair_set],
        db=IN_SILICO_DIR / "indexes/{genotype}/genome.fasta",
        index_done=IN_SILICO_DIR / "indexes/{genotype}/index.done",

    output:
        report=IN_SILICO_DIR / "specificity/{genotype}/{pair_set}",
        mfe_log=IN_SILICO_DIR / "specificity/{genotype}/{pair_set}.mfe.log",
        spec=IN_SILICO_DIR / "specificity/{genotype}/{pair_set}.spec.tsv",

    params:
        out_prefix=lambda wildcards: (
            IN_SILICO_DIR
            / "specificity"
            / wildcards.genotype
            / wildcards.pair_set
        ),
        min_tm=lambda wildcards: (
            config.get("kasp", {})
            .get("mfeprimer", {})
            .get("specificity", {})
            .get("min_tm", 40)
        ),
        extra=lambda wildcards: (
            config.get("advanced", {})
            .get("mfeprimer", {})
            .get("specificity_extra_options", [])
        ),

    threads: 1

    log:
        stdout=LOG_DIR / "run_mfeprimer_specificity/{genotype}.{pair_set}.stdout",
        stderr=LOG_DIR / "run_mfeprimer_specificity/{genotype}.{pair_set}.stderr",

    benchmark:
        BENCHMARK_DIR / "run_mfeprimer_specificity/{genotype}.{pair_set}.tsv"

    shell:
        r"""
        mkdir -p "$(dirname "{output.spec}")" "$(dirname "{log.stdout}")"

        if [[ ! -s "{input.primers}" ]]; then
            touch "{output.report}" "{output.mfe_log}"
            printf '#1-based coordinate, should be idential to blat output if align primer sequences to genome\n#name\tchrom\tampStart\tampEnd\tampGC\tampSize\tfpName\tfpStart\tfpEnd\tfpSeq\tfpTm\tfpGC\tfpDg\trpName\trpEnd\trpStart\trpSeq\trpTm\trpGC\trpDg\tnote\n' \
                > "{output.spec}"
        else
            "{workflow.basedir}/dependencies/bin/mfeprimer-4.2.3-linux-amd64" spec \
                --tm "{params.min_tm}" \
                --cpu {threads} \
                -i "{input.primers}" \
                -d "{input.db}" \
                -o "{output.report}" \
                {params.extra:q} \
                1> "{log.stdout}" \
                2> "{log.stderr}"

            if [[ ! -f "{output.spec}" ]]; then
                printf '#1-based coordinate, should be idential to blat output if align primer sequences to genome\n#name\tchrom\tampStart\tampEnd\tampGC\tampSize\tfpName\tfpStart\tfpEnd\tfpSeq\tfpTm\tfpGC\tfpDg\trpName\trpEnd\trpStart\trpSeq\trpTm\trpGC\trpDg\tnote\n' \
                    > "{output.spec}"
            fi
        fi
        """

rule run_mfeprimer_dimer:
    input:
        primers=POLYMARKER_SUMMARY_DIR / "primers_with_tails.fasta",

    output:
        report=IN_SILICO_DIR / "dimers.tsv",

    params:
        max_dg=lambda wildcards: (
            config.get("kasp", {})
            .get("mfeprimer", {})
            .get("dimer", {})
            .get("max_dg", -3.5)
        ),
        extra=lambda wildcards: (
            config.get("advanced", {})
            .get("mfeprimer", {})
            .get("dimer_extra_options", [])
        ),

    threads: 1

    log:
        stdout=LOG_DIR / "run_mfeprimer_dimer.stdout",
        stderr=LOG_DIR / "run_mfeprimer_dimer.stderr",

    benchmark:
        BENCHMARK_DIR / "run_mfeprimer_dimer.tsv"

    shell:
        r"""
        mkdir -p "$(dirname "{output.report}")" "$(dirname "{log.stdout}")"

        if [[ ! -s "{input.primers}" ]]; then
            touch "{output.report}"
        else
            "{workflow.basedir}/dependencies/bin/mfeprimer-4.2.3-linux-amd64" dimer \
                --primer \
                --dg "{params.max_dg}" \
                --cpu {threads} \
                -i "{input.primers}" \
                -o "{output.report}" \
                {params.extra:q} \
                1> "{log.stdout}" \
                2> "{log.stderr}"
        fi
        """


rule run_mfeprimer_hairpin:
    input:
        primers=POLYMARKER_SUMMARY_DIR / "primers_with_tails.fasta",

    output:
        report=IN_SILICO_DIR / "hairpins.tsv",

    params:
        extra=lambda wildcards: (
            config.get("advanced", {})
            .get("mfeprimer", {})
            .get("hairpin_extra_options", [])
        ),

    threads: 1

    log:
        stdout=LOG_DIR / "run_mfeprimer_hairpin.stdout",
        stderr=LOG_DIR / "run_mfeprimer_hairpin.stderr",

    benchmark:
        BENCHMARK_DIR / "run_mfeprimer_hairpin.tsv"

    shell:
        r"""
        mkdir -p "$(dirname "{output.report}")" "$(dirname "{log.stdout}")"

        if [[ ! -s "{input.primers}" ]]; then
            : > "{output.report}"
        else
            "{workflow.basedir}/dependencies/bin/mfeprimer-4.2.3-linux-amd64" hairpin \
                --cpu {threads} \
                -i "{input.primers}" \
                -o "{output.report}" \
                {params.extra:q} \
                1> "{log.stdout}" \
                2> "{log.stderr}"
        fi
        """


rule summarize_in_silico_validation:
    input:
        snp_status=POLYMARKER_SUMMARY_DIR / "polymarker_snp_status.tsv",
        snp_status_by_genotype=(
            POLYMARKER_SUMMARY_DIR
            / "polymarker_snp_status_by_genotype.tsv"
        ),
        assays=POLYMARKER_SUMMARY_DIR / "polymarker_assays.tsv",
        snp_positions=SNP_POS_LONG_TSV,
        specificity=expand(
            str(
                IN_SILICO_DIR
                / "specificity/{genotype}/{pair_set}.spec.tsv"
            ),
            genotype=KASP_GENOTYPE_NAMES,
            pair_set=["canonical", "noncanonical"],
        ),
        dimers=IN_SILICO_DIR / "dimers.tsv",
        hairpins=IN_SILICO_DIR / "hairpins.tsv",

    output:
        assay_status=VALIDATION_DIR / "assay_status.tsv",
        assay_status_by_genotype=(
            VALIDATION_DIR / "assay_status_by_genotype.tsv"
        ),
        snp_status=VALIDATION_DIR / "snp_status.tsv",
        validated_assays=VALIDATION_DIR / "validated_assays.tsv",

    log:
        stdout=LOG_DIR / "summarize_in_silico_validation.stdout",
        stderr=LOG_DIR / "summarize_in_silico_validation.stderr",

    benchmark:
        BENCHMARK_DIR / "summarize_in_silico_validation.tsv"

    shell:
        r"""
        mkdir -p "{VALIDATION_DIR}" "$(dirname "{log.stdout}")"

        python3 "{SCRIPTS_DIR}/summarize_in_silico_validation.py" \
            --polymarker-summary "{POLYMARKER_SUMMARY_DIR}" \
            --snp-positions "{input.snp_positions}" \
            --in-silico-dir "{IN_SILICO_DIR}" \
            --output-dir "{VALIDATION_DIR}" \
            1> "{log.stdout}" \
            2> "{log.stderr}"
        """
