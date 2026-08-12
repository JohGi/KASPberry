import polars as pl


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
