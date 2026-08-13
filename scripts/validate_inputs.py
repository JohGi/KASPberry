from pathlib import Path

import polars as pl
from snakemake.utils import validate

GENOTYPE_COLUMNS = {
    "genotype",
    "group",
    "region_fasta",
    "genome_fasta",
    "source_seq",
    "region_start",
}

ANNOTATION_COLUMNS = {
    "genotype",
    "track",
    "gff",
}


def check_genotypes_table(df: pl.DataFrame) -> None:
    """Validate constraints involving the genotype table as a whole."""

    if df.height == 0:
        raise ValueError("Genotypes table must not be empty")

    missing = GENOTYPE_COLUMNS - set(df.columns)
    if missing:
        raise ValueError(
            "Missing required columns in genotypes.tsv: "
            + ", ".join(sorted(missing))
        )

    duplicates = (
        df.group_by("genotype")
        .len()
        .filter(pl.col("len") > 1)
        .get_column("genotype")
        .to_list()
    )

    if duplicates:
        raise ValueError(
            "Duplicate genotype names in genotypes.tsv: "
            + ", ".join(sorted(duplicates))
        )

    groups = (
        df.get_column("group")
        .drop_nulls()
        .unique()
        .to_list()
    )

    if len(groups) not in {0, 2}:
        raise ValueError(
            "Column 'group' must contain either no groups or exactly "
            f"two distinct groups. Found: {sorted(groups)}"
        )

    has_genome = pl.col("genome_fasta").is_not_null()
    has_source = pl.col("source_seq").is_not_null()
    has_start = pl.col("region_start").is_not_null()
    incomplete_metadata = df.filter(
        (has_genome & ~(has_source & has_start))
        | (has_source ^ has_start)
    )
    if incomplete_metadata.height:
        names = incomplete_metadata.get_column("genotype").to_list()
        raise ValueError(
            "genome_fasta, source_seq, and region_start must be provided "
            "together in genotypes.tsv for: "
            + ", ".join(names)
        )


def check_annotations_table(
    df: pl.DataFrame,
    genotype_names: set[str],
) -> None:
    """Validate constraints involving annotations and genotype names."""

    missing = ANNOTATION_COLUMNS - set(df.columns)
    if missing:
        raise ValueError(
            "Missing required columns in annotations.tsv: "
            + ", ".join(sorted(missing))
        )

    unknown = (
        df.filter(~pl.col("genotype").is_in(genotype_names))
        .get_column("genotype")
        .unique()
        .to_list()
    )

    if unknown:
        raise ValueError(
            "Unknown genotypes in annotations.tsv: "
            + ", ".join(sorted(unknown))
        )


def check_chromosomes_table(
    chromosomes: pl.DataFrame,
    genotypes: pl.DataFrame,
    polymarker_genomes: int,
) -> None:
    """Check semantic consistency of chromosomes.tsv."""

    # (genotype, seq_id) must be unique.
    duplicates = (
        chromosomes
        .group_by(["genotype", "seq_id"])
        .len()
        .filter(pl.col("len") > 1)
    )

    if duplicates.height:
        duplicate_pairs = [
            f"{row['genotype']} / {row['seq_id']}"
            for row in duplicates.iter_rows(named=True)
        ]
        raise ValueError(
            "Duplicate (genotype, seq_id) entries in chromosomes.tsv: "
            + ", ".join(sorted(duplicate_pairs))
        )

    # All genotypes in chromosomes.tsv must exist in genotypes.tsv.
    unknown_genotypes = (
        set(chromosomes.get_column("genotype").to_list())
        - set(genotypes.get_column("genotype").to_list())
    )

    if unknown_genotypes:
        raise ValueError(
            "Unknown genotypes in chromosomes.tsv: "
            + ", ".join(sorted(unknown_genotypes))
        )

    # Only genotypes with complete whole-genome metadata are KASP QC genomes.
    qc_genotypes = genotypes.filter(
        pl.col("genome_fasta").is_not_null()
        & pl.col("source_seq").is_not_null()
        & pl.col("region_start").is_not_null()
    )

    for genotype_row in qc_genotypes.iter_rows(named=True):
        genotype = genotype_row["genotype"]
        source_seq = genotype_row["source_seq"]

        genotype_chromosomes = chromosomes.filter(
            pl.col("genotype") == genotype
        )

        source_rows = genotype_chromosomes.filter(
            pl.col("seq_id") == source_seq
        )

        if source_rows.height == 0:
            raise ValueError(
                f"source_seq '{source_seq}' for genotype '{genotype}' "
                "is missing from chromosomes.tsv."
            )

        target_group = source_rows.get_column(
            "homoeologous_group"
        )[0]

        group_rows = genotype_chromosomes.filter(
            pl.col("homoeologous_group") == target_group
        )

        n_subgenomes = (
            group_rows
            .get_column("subgenome")
            .n_unique()
        )

        if n_subgenomes != polymarker_genomes:
            raise ValueError(
                f"Homoeologous group '{target_group}' for genotype "
                f"'{genotype}' contains {n_subgenomes} distinct subgenomes, "
                f"but kasp.polymarker_genomes={polymarker_genomes}."
            )


def validate_inputs(
    config: dict,
    mode: str,
    schema_dir: str | Path,
) -> None:
    """Validate KASPberry input tables for the requested workflow mode."""

    if mode not in {"snps", "kasp"}:
        raise ValueError(
            f"Unknown workflow mode '{mode}'. Expected 'snps' or 'kasp'."
        )

    schema_dir = Path(schema_dir)

    # ------------------------------------------------------------------
    # Configuration
    # ------------------------------------------------------------------

    validate(
        config,
        schema_dir / "config.schema.yaml",
        set_default=False,
    )

    # ------------------------------------------------------------------
    # genotypes.tsv
    # Required by both workflows.
    # ------------------------------------------------------------------

    genotypes_path = Path(config["inputs"]["genotypes"])

    if not genotypes_path.is_file():
        raise ValueError(
            f"genotypes.tsv not found: {genotypes_path}"
        )

    genotypes_df = pl.read_csv(
        genotypes_path,
        separator="\t",
        null_values="",
    )

    validate(
        genotypes_df,
        schema_dir / "genotypes.schema.yaml",
        set_default=False,
    )

    check_genotypes_table(genotypes_df)

    genotype_names = set(
        genotypes_df.get_column("genotype").to_list()
    )

    # ------------------------------------------------------------------
    # annotations.tsv
    # Optional for both workflows.
    # ------------------------------------------------------------------

    annotations_path = config.get("inputs", {}).get("annotations")

    if annotations_path:
        annotations_path = Path(annotations_path)

        if not annotations_path.is_file():
            raise ValueError(
                f"annotations.tsv not found: {annotations_path}"
            )

        annotations_df = pl.read_csv(
            annotations_path,
            separator="\t",
            null_values="",
        )

        validate(
            annotations_df,
            schema_dir / "annotations.schema.yaml",
            set_default=False,
        )

        check_annotations_table(
            annotations_df,
            genotype_names=genotype_names,
        )

    # Nothing below is required for SNP discovery alone.
    if mode == "snps":
        return

    # ------------------------------------------------------------------
    # KASP-specific inputs
    # ------------------------------------------------------------------

    kasp_config = config.get("kasp")

    if kasp_config is None:
        raise ValueError(
            "Section 'kasp' is required when running the KASP workflow."
        )

    if "polymarker_genomes" not in kasp_config:
        raise ValueError(
            "kasp.polymarker_genomes is required when running "
            "the KASP workflow."
        )

    polymarker_genomes = kasp_config["polymarker_genomes"]

    # At least one genotype must provide a whole genome for KASP QC.
    qc_genotypes_df = genotypes_df.filter(
        pl.col("genome_fasta").is_not_null()
    )

    if qc_genotypes_df.height == 0:
        raise ValueError(
            "The KASP workflow requires at least one genotype with "
            "genome_fasta, source_seq, and region_start."
        )

    # ------------------------------------------------------------------
    # chromosomes.tsv
    #
    # Required for PolyMarker when several subgenomes must be distinguished.
    # Optional when polymarker_genomes == 1.
    # ------------------------------------------------------------------

    chromosomes_path = config.get("inputs", {}).get("chromosomes")

    if polymarker_genomes > 1 and not chromosomes_path:
        raise ValueError(
            "inputs.chromosomes is required when "
            "kasp.polymarker_genomes > 1."
        )

    if chromosomes_path:
        chromosomes_path = Path(chromosomes_path)

        if not chromosomes_path.is_file():
            raise ValueError(
                f"chromosomes.tsv not found: {chromosomes_path}"
            )

        chromosomes_df = pl.read_csv(
            chromosomes_path,
            separator="\t",
            null_values="",
            infer_schema=False,
        )

        validate(
            chromosomes_df,
            schema_dir / "chromosomes.schema.yaml",
            set_default=False,
        )

        check_chromosomes_table(
            chromosomes_df,
            genotypes=genotypes_df,
            polymarker_genomes=polymarker_genomes,
        )
