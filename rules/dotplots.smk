from itertools import combinations


def resolve_dotplot_pairs(
    sample_names: list[str],
    config: dict,
) -> list[tuple[str, str]]:
    """Resolve pairwise dotplot comparisons from config."""
    pivot = str(config["viewer"]["dotplot_reference"] or "").strip()

    if not pivot:
        return list(combinations(sample_names, 2))

    if pivot not in sample_names:
        raise ValueError(
            f"Unknown viewer.dotplot_reference: {pivot!r}. "
            f"Expected one of: {sample_names}"
        )

    return [(pivot, sample) for sample in sample_names if sample != pivot]


def build_pair_id(sample_a: str, sample_b: str) -> str:
    """Build a stable pair identifier."""
    return f"{sample_a}__vs__{sample_b}"


def split_pair_id(pair_id: str) -> tuple[str, str]:
    """Decode a pair identifier."""
    parts = pair_id.split("__vs__")
    if len(parts) != 2 or not parts[0] or not parts[1]:
        raise ValueError(f"Invalid pair_id: {pair_id!r}")
    return parts[0], parts[1]


def get_pair_sample_a(wildcards) -> str:
    sample_a, _sample_b = split_pair_id(wildcards.pair_id)
    return sample_a


def get_pair_sample_b(wildcards) -> str:
    _sample_a, sample_b = split_pair_id(wildcards.pair_id)
    return sample_b


DOTPLOT_PAIRS = resolve_dotplot_pairs(GENOTYPE_NAMES, config)
DOTPLOT_PAIR_IDS = [build_pair_id(sample_a, sample_b) for sample_a, sample_b in DOTPLOT_PAIRS]
DOTPLOT_PAFS = expand(DOTPLOT_PAF_DIR / "{pair_id}.paf", pair_id=DOTPLOT_PAIR_IDS)
DOTPLOT_FORMATTED = expand(DOTPLOT_FORMATTED_DIR / "{pair_id}.tsv", pair_id=DOTPLOT_PAIR_IDS)
DOTPLOT_SVGS = expand(DOTPLOT_SVG_DIR / "{pair_id}.svg", pair_id=DOTPLOT_PAIR_IDS)
DOTPLOT_ONLY_SVGS = expand(
    DOTPLOT_ONLY_SVG_DIR / "{pair_id}.dotplot_only.svg",
    pair_id=DOTPLOT_PAIR_IDS,
)


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
            config["viewer"]["dotplot_reference"] or ""
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
