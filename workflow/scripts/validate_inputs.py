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


def resolve_input_path(
    path: str | Path,
    working_directory: Path,
) -> Path:
    """Resolve a user input path as Snakemake does from its workdir."""
    input_path = Path(path).expanduser()
    if input_path.is_absolute():
        return input_path.resolve()
    return (working_directory / input_path).resolve()


def _iter_fasta_records(path: Path):
    """Yield ``(record_id, sequence)`` pairs from a FASTA file."""
    try:
        handle = path.open("r", encoding="utf-8")
    except OSError as error:
        raise ValueError(f"FASTA file is not readable: {path}: {error}") from error

    current_id = None
    sequence_parts = []
    try:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    yield current_id, "".join(sequence_parts)
                header = line[1:].strip()
                if not header:
                    raise ValueError(f"FASTA record has an empty identifier: {path}")
                current_id = header.split()[0]
                sequence_parts = []
            else:
                if current_id is None:
                    raise ValueError(f"Invalid FASTA file (sequence before header): {path}")
                sequence_parts.append("".join(line.split()))
        if current_id is not None:
            yield current_id, "".join(sequence_parts)
    finally:
        handle.close()


def _read_single_region_fasta(path: Path) -> str:
    """Validate a regional FASTA and return its sequence."""
    records = list(_iter_fasta_records(path))
    if len(records) != 1 or not records[0][1]:
        raise ValueError(
            f"region_fasta must contain exactly one non-empty FASTA record: {path}"
        )
    return records[0][1]


def _inspect_genome_fasta(
    path: Path,
    target_id: str | None = None,
    segment_start: int | None = None,
    segment_length: int = 0,
) -> tuple[set[str], int | None, str]:
    """Stream a genome FASTA, retaining only a requested target segment."""
    ids: set[str] = set()
    current_id = None
    current_length = 0
    target_length = None
    segment_parts: list[str] = []
    segment_end = (
        segment_start + segment_length
        if segment_start is not None
        else None
    )

    try:
        handle = path.open("r", encoding="utf-8")
    except OSError as error:
        raise ValueError(f"genome_fasta is not readable: {path}: {error}") from error

    try:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None and current_id == target_id:
                    target_length = current_length
                header = line[1:].strip()
                if not header:
                    raise ValueError(f"FASTA record has an empty identifier: {path}")
                current_id = header.split()[0]
                if current_id in ids:
                    raise ValueError(f"Duplicate FASTA record identifier '{current_id}' in {path}")
                ids.add(current_id)
                current_length = 0
                continue

            if current_id is None:
                raise ValueError(f"Invalid FASTA file (sequence before header): {path}")
            sequence = "".join(line.split())
            if current_id == target_id and segment_start is not None:
                overlap_start = max(segment_start, current_length)
                overlap_end = min(segment_end, current_length + len(sequence))
                if overlap_start < overlap_end:
                    start = overlap_start - current_length
                    end = overlap_end - current_length
                    segment_parts.append(sequence[start:end])
            current_length += len(sequence)

        if current_id is not None and current_id == target_id:
            target_length = current_length
    finally:
        handle.close()

    return ids, target_length, "".join(segment_parts)


def _validate_shared_files(
    genotypes: pl.DataFrame,
    annotations: pl.DataFrame | None,
    config: dict,
    genotype_names: set[str],
    working_directory: Path,
) -> dict[str, str]:
    """Validate files and references shared by SNP and KASP workflows."""

    region_sequences: dict[str, str] = {}

    for row in genotypes.iter_rows(named=True):
        path = resolve_input_path(
            row["region_fasta"],
            working_directory,
        )

        if not path.is_file():
            raise ValueError(f"region_fasta not found: {path}")

        region_sequences[row["genotype"]] = _read_single_region_fasta(path)

    if annotations is not None:
        for gff in annotations.get_column("gff").unique().to_list():
            gff_path = resolve_input_path(
                gff,
                working_directory,
            )
            if not gff_path.is_file():
                raise ValueError(f"GFF file not found: {gff}")

    library = (
        config
        .get("snps", {})
        .get("repeat_masking", {})
        .get("library")
    )

    if library and not resolve_input_path(
        library,
        working_directory,
    ).is_file():
        raise ValueError(
            f"Repeat-masking library not found: {library}"
        )

    dotplot_reference = (
        config
        .get("viewer", {})
        .get("dotplot_reference")
    )

    if dotplot_reference and dotplot_reference not in genotype_names:
        raise ValueError(
            "Unknown viewer.dotplot_reference: "
            f"{dotplot_reference!r}. "
            "Expected a genotype from genotypes.tsv."
        )

    return region_sequences


def _read_gff_seqids(path: Path) -> set[str]:
    """Return seqids from non-comment GFF records."""
    seqids: set[str] = set()

    try:
        handle = path.open("r", encoding="utf-8")
    except OSError as error:
        raise ValueError(f"GFF file is not readable: {path}: {error}") from error

    try:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.rstrip("\r\n")
            if stripped.strip() == "##FASTA":
                break
            if not stripped.strip() or stripped.lstrip().startswith("#"):
                continue

            fields = stripped.split("\t")
            if len(fields) != 9:
                raise ValueError(
                    f"Invalid GFF line in {path} at line {line_number}: "
                    f"expected 9 columns, got {len(fields)}."
                )
            seqids.add(fields[0])
    finally:
        handle.close()

    return seqids


def check_annotation_source_sequences(
    genotypes: pl.DataFrame,
    annotations: pl.DataFrame,
    working_directory: Path,
) -> None:
    """Require GFF-backed genotypes to declare a matching source sequence."""
    source_seq_by_genotype = dict(
        genotypes.select("genotype", "source_seq").iter_rows()
    )
    gff_seqids_cache: dict[Path, set[str]] = {}

    for row in annotations.iter_rows(named=True):
        genotype = row["genotype"]
        source_seq = source_seq_by_genotype[genotype]
        if source_seq is None:
            raise ValueError(
                f"Genotype '{genotype}' has a GFF configured but is missing "
                "source_seq and region_start."
            )

        gff_path = resolve_input_path(
            row["gff"],
            working_directory,
        )
        if gff_path not in gff_seqids_cache:
            gff_seqids_cache[gff_path] = _read_gff_seqids(gff_path)

        if source_seq not in gff_seqids_cache[gff_path]:
            raise ValueError(
                f"GFF '{gff_path}' for genotype '{genotype}' does not contain "
                f"source_seq '{source_seq}'."
            )


def _validate_kasp_genomes(
    genotypes: pl.DataFrame,
    region_sequences: dict[str, str],
    working_directory: Path,
) -> dict[str, set[str]]:
    """Validate KASP genome coordinates and return FASTA IDs by genotype."""
    genome_ids: dict[str, set[str]] = {}
    for row in genotypes.iter_rows(named=True):
        genome_path = row["genome_fasta"]
        if genome_path is None:
            continue
        genotype = row["genotype"]
        start = int(row["region_start"])
        region = region_sequences[genotype]
        ids, source_length, prefix = _inspect_genome_fasta(
            resolve_input_path(genome_path, working_directory),
            target_id=row["source_seq"],
            segment_start=start - 1,
            segment_length=min(500, len(region)),
        )
        genome_ids[genotype] = ids
        if source_length is None:
            raise ValueError(
                f"source_seq '{row['source_seq']}' for genotype '{genotype}' "
                f"was not found in genome_fasta: {genome_path}"
            )
        region_end = start + len(region) - 1
        if region_end > source_length:
            raise ValueError(
                f"Region for genotype '{genotype}' extends to {region_end}, "
                f"beyond source_seq '{row['source_seq']}' length {source_length}."
            )
        expected = region[: min(500, len(region))].upper()
        if prefix.upper() != expected:
            raise ValueError(
                f"The first {len(expected)} bases of region_fasta for genotype "
                f"'{genotype}' do not match source_seq '{row['source_seq']}' "
                f"at 1-based region_start {start}."
            )
    return genome_ids


def check_extra_options(
    options: list[str],
    forbidden: set[str],
    config_key: str,
) -> None:
    """Reject command-line options managed directly by KASPberry."""
    for option in options:
        flag = option.split("=", 1)[0]

        if flag in forbidden:
            raise ValueError(
                f"{config_key} must not redefine KASPberry-managed "
                f"option '{flag}'."
            )


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

    incomplete_source_coordinates = df.filter(has_source != has_start)
    if incomplete_source_coordinates.height:
        names = incomplete_source_coordinates.get_column("genotype").to_list()
        raise ValueError(
            "source_seq and region_start must be provided together in "
            "genotypes.tsv for: "
            + ", ".join(names)
        )

    incomplete_genome_coordinates = df.filter(
        has_genome & ~(has_source & has_start)
    )
    if incomplete_genome_coordinates.height:
        names = incomplete_genome_coordinates.get_column("genotype").to_list()
        raise ValueError(
            "genome_fasta requires source_seq and region_start in "
            "genotypes.tsv for: "
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

    duplicates = (
        df.group_by(["genotype", "track"])
        .len()
        .filter(pl.col("len") > 1)
    )

    if duplicates.height:
        pairs = [
            f"{row['genotype']} / {row['track']}"
            for row in duplicates.iter_rows(named=True)
        ]
        raise ValueError(
            "Duplicate (genotype, track) entries in annotations.tsv: "
            + ", ".join(sorted(pairs))
        )


def check_chromosomes_table(
    chromosomes: pl.DataFrame,
    genotypes: pl.DataFrame,
    polymarker_genomes: int,
    genome_ids_by_genotype: dict[str, set[str]] | None = None,
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

    if genome_ids_by_genotype is not None:
        no_genome = set(chromosomes.get_column("genotype").to_list()) - set(
            genome_ids_by_genotype
        )
        if no_genome:
            raise ValueError(
                "Chromosome rows require a genome_fasta for genotype(s): "
                + ", ".join(sorted(no_genome))
            )
        missing_seq_ids = []
        for row in chromosomes.iter_rows(named=True):
            if row["seq_id"] not in genome_ids_by_genotype[row["genotype"]]:
                missing_seq_ids.append(
                    f"{row['genotype']} / {row['seq_id']}"
                )
        if missing_seq_ids:
            raise ValueError(
                "Unknown seq_id in chromosomes.tsv for genome_fasta: "
                + ", ".join(sorted(missing_seq_ids))
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
    working_directory: str | Path,
) -> None:
    """Validate KASPberry input tables for the requested workflow mode."""

    if mode not in {"snps", "kasp"}:
        raise ValueError(
            f"Unknown workflow mode '{mode}'. Expected 'snps' or 'kasp'."
        )

    schema_dir = Path(schema_dir)
    working_directory = Path(working_directory).expanduser().resolve()

    # ------------------------------------------------------------------
    # Configuration
    # ------------------------------------------------------------------

    validate(
        config,
        schema_dir / "config.schema.yaml",
    )
    
    if mode == "kasp" and config["snps"]["min_snp_flank"] < 50:
        raise ValueError(
            "snps.min_snp_flank must be at least 50 for the KASP workflow."
        )

    # ------------------------------------------------------------------
    # genotypes.tsv
    # Required by both workflows.
    # ------------------------------------------------------------------

    genotypes_path = resolve_input_path(
        config["inputs"]["genotypes"],
        working_directory,
    )

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

    annotations_df = None
    if annotations_path:
        annotations_path = resolve_input_path(
            annotations_path,
            working_directory,
        )

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

    region_sequences = _validate_shared_files(
        genotypes_df,
        annotations_df,
        config,
        genotype_names,
        working_directory,
    )

    if annotations_df is not None:
        check_annotation_source_sequences(
            genotypes_df,
            annotations_df,
            working_directory,
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

    genome_ids_by_genotype = _validate_kasp_genomes(
        genotypes_df,
        region_sequences,
        working_directory,
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
        chromosomes_path = resolve_input_path(
            chromosomes_path,
            working_directory,
        )

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
            genome_ids_by_genotype=genome_ids_by_genotype,
        )

    # ------------------------------------------------------------------
    # MFE-primer advanced options
    # ------------------------------------------------------------------

    mfeprimer_extra = (
        config.get("advanced", {})
        .get("mfeprimer", {})
    )

    check_extra_options(
        mfeprimer_extra.get("specificity_extra_options", []),
        {
            "-i", "--in",
            "-o", "--out",
            "-d", "--db",
            "-c", "--cpu",
            "-t", "--tm",
            "-j", "--json",
            "-g", "--gz",
            "-f", "--fasta",
            "--PE",
            "--SE",
            "--cutprimer",
            "--virus",
        },
        "advanced.mfeprimer.specificity_extra_options",
    )

    check_extra_options(
        mfeprimer_extra.get("dimer_extra_options", []),
        {
            "-i", "--in",
            "-o", "--out",
            "-c", "--cpu",
            "-p", "--primer",
            "-d", "--dg",
            "-j", "--json",
        },
        "advanced.mfeprimer.dimer_extra_options",
    )

    check_extra_options(
        mfeprimer_extra.get("hairpin_extra_options", []),
        {
            "-i", "--in",
            "-o", "--out",
            "-c", "--cpu",
            "-j", "--json",
        },
        "advanced.mfeprimer.hairpin_extra_options",
    )
