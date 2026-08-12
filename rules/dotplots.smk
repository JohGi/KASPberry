rule run_pairwise_minimap2:
    input:
        fasta_a=lambda wildcards: CLEAN_FASTA_DIR / f"{get_pair_sample_a(wildcards)}.fasta",
        fasta_b=lambda wildcards: CLEAN_FASTA_DIR / f"{get_pair_sample_b(wildcards)}.fasta"
    output:
        DOTPLOT_PAF_DIR / "{pair_id}.paf"
    benchmark:
        BENCHMARK_DIR / "run_pairwise_minimap2" / "{pair_id}.tsv"
    log:
        LOG_DIR / "run_pairwise_minimap2" / "{pair_id}.stderr"
    threads: 1
    shell:
        r"""
        mkdir -p "{DOTPLOT_PAF_DIR}" "$(dirname "{log}")"
        minimap2 \
            -x asm5 \
            -c \
            -t {threads} \
            "{input.fasta_a}" \
            "{input.fasta_b}" \
            > "{output}" \
            2> "{log}"
        """

rule format_pairwise_paf_for_blastn2dotplots:
    input:
        DOTPLOT_PAF_DIR / "{pair_id}.paf"
    output:
        DOTPLOT_FORMATTED_DIR / "{pair_id}.tsv"
    benchmark:
        BENCHMARK_DIR / "format_pairwise_paf_for_blastn2dotplots" / "{pair_id}.tsv"
    log:
        LOG_DIR / "format_pairwise_paf_for_blastn2dotplots" / "{pair_id}.stderr"
    params:
        converter=SCRIPTS_DIR / "blastn2dotplots/utilities/paf2blastn-fmt6.pl"
    shell:
        r"""
        mkdir -p "{DOTPLOT_FORMATTED_DIR}" "$(dirname "{log}")"
        perl "{params.converter}" "{input}" > "{output}" 2> "{log}"
        """

rule run_pairwise_blastn2dotplots:
    input:
        formatted=DOTPLOT_FORMATTED_DIR / "{pair_id}.tsv",
        fasta_a=lambda wildcards: CLEAN_FASTA_DIR / f"{get_pair_sample_a(wildcards)}.fasta",
        fasta_b=lambda wildcards: CLEAN_FASTA_DIR / f"{get_pair_sample_b(wildcards)}.fasta"
    output:
        standard=DOTPLOT_PDF_DIR / "{pair_id}.pdf",
        dotplot_only=DOTPLOT_PDF_DIR / "{pair_id}.dotplot_only.pdf",
    benchmark:
        BENCHMARK_DIR / "run_pairwise_blastn2dotplots" / "{pair_id}.tsv"
    log:
        stdout=LOG_DIR / "run_pairwise_blastn2dotplots" / "{pair_id}.stdout",
        stderr=LOG_DIR / "run_pairwise_blastn2dotplots" / "{pair_id}.stderr"
    params:
        db_name=get_pair_sample_a,
        query_name=get_pair_sample_b,
        out_prefix=lambda wildcards: DOTPLOT_PDF_DIR / f"{wildcards.pair_id}",
    shell:
        r"""
        mkdir -p "{DOTPLOT_PDF_DIR}" "$(dirname "{log.stdout}")"

        pixi run -e dotplot bash "{SCRIPTS_DIR}/run_pairwise_blastn2dotplots.sh" \
            --blastn-tsv "{input.formatted}" \
            --db-name "{params.db_name}" \
            --query-name "{params.query_name}" \
            --out-prefix "{params.out_prefix}" \
            1> "{log.stdout}" \
            2> "{log.stderr}"
        """

rule convert_pairwise_dotplot_pdf_to_svg:
    input:
        DOTPLOT_PDF_DIR / "{pair_id}.pdf"
    output:
        svg=DOTPLOT_SVG_DIR / "{pair_id}.svg"
    benchmark:
        BENCHMARK_DIR / "convert_pairwise_dotplot_pdf_to_svg" / "{pair_id}.tsv"
    log:
        LOG_DIR / "convert_pairwise_dotplot_pdf_to_svg" / "{pair_id}.stderr"
    shell:
        r"""
        mkdir -p "{DOTPLOT_SVG_DIR}" "$(dirname "{log}")"

        pdf2svg "{input}" "{output.svg}" 2> "{log}"
        """

rule convert_pairwise_dotplot_only_pdf_to_svg:
    input:
        DOTPLOT_PDF_DIR / "{pair_id}.dotplot_only.pdf"
    output:
        svg=DOTPLOT_ONLY_SVG_DIR / "{pair_id}.dotplot_only.svg"
    benchmark:
        BENCHMARK_DIR / "convert_pairwise_dotplot_only_pdf_to_svg" / "{pair_id}.tsv"
    log:
        LOG_DIR / "convert_pairwise_dotplot_only_pdf_to_svg" / "{pair_id}.stderr"
    shell:
        r"""
        mkdir -p "{DOTPLOT_ONLY_SVG_DIR}" "$(dirname "{log}")"

        pdf2svg "{input}" "{output.svg}" 2> "{log}"
        """

rule build_dotplot_gallery_html:
    input:
        svgs=DOTPLOT_SVGS
    output:
        html=DOTPLOT_GALLERY_HTML
    benchmark:
        BENCHMARK_DIR / "build_dotplot_gallery_html" / "dotplots_gallery.tsv"
    log:
        stdout=LOG_DIR / "build_dotplot_gallery_html" / "dotplots_gallery.stdout",
        stderr=LOG_DIR / "build_dotplot_gallery_html" / "dotplots_gallery.stderr"
    params:
        pivot=lambda wildcards: str(
            config.get("viewer", {}).get("dotplot_reference") or ""
        ).strip()
    shell:
        r"""
        mkdir -p "{DOTPLOT_COMBINED_DIR}" "$(dirname "{log.stdout}")"
        python3 "{SCRIPTS_DIR}/build_dotplot_gallery_html.py" \
            --samples "{GENOTYPES_TSV}" \
            --svg-dir "{DOTPLOT_SVG_DIR}" \
            --output "{output.html}" \
            --pivot "{params.pivot}" \
            --title "{PROJECT_TITLE}" \
            1> "{log.stdout}" \
            2> "{log.stderr}"
        """
