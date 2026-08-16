#!/usr/bin/env python3
"""Build a synthetic homologous test dataset from block FASTA files."""

import argparse
import random
import re
from pathlib import Path

from attrs import define, field
from attrs.validators import instance_of

VALID_FASTA_EXTENSIONS = {".fa", ".fasta", ".fna"}
DNA_ALPHABET = "ACGT"
FASTA_LINE_WIDTH = 80

GENOTYPE_NAMES = ["Karur", "Soldur", "Belalur", "Lloyd", "Dic2"]

FASTA_HEADER_ALIASES = {
    "Dic2": "Dic2_mlk",
}

# Synthetic chromosome metadata used by chromosomes.tsv.
#
# Each generated "whole-genome" FASTA contains the QTL-bearing pseudochromosome
# (subgenome A) plus one synthetic homoeologue (subgenome B). The B sequence is
# derived from A by deterministic substitutions only, so both records have the
# same length and coordinate system while remaining highly similar.
SYNTHETIC_HOMOEOLOGOUS_GROUP = "7"
SYNTHETIC_TARGET_SUBGENOME = "A"
SYNTHETIC_HOMOEOLOGUE_SUBGENOME = "B"
SYNTHETIC_HOMOEOLOGUE_SUFFIX = "_B"
SYNTHETIC_HOMOEOLOGUE_MUTATION_INTERVAL = 50
SYNTHETIC_HOMOEOLOGUE_MUTATION_OFFSET = 25

BLOCK_ORDER = [
    "conserved_4_not_5",
    "too_short",
    "mixed_strand",
    "two_good_snps_dic2",
    "two_good_snps_soldur",
    "too_close_to_masked",
    "too_close_snps",
    "too_close_to_indel",
    "too_close_to_start",
]


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Build a synthetic QTL dataset from block FASTA files. "
            "The script writes one QTL FASTA per genotype, one pseudo-chromosome "
            "FASTA per genotype, a genotypes.tsv file, a chromosomes.tsv file, "
            "and summary TSV files."
        )
    )
    parser.add_argument(
        "--input-dir",
        required=True,
        help="Directory containing one FASTA file per block.",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Output directory.",
    )
    parser.add_argument(
        "--spacer-length",
        default="300",
        help=(
            "Spacer length specification. Use a single integer for fixed length "
            "(e.g. 300) or a range MIN-MAX for random length (e.g. 250-350)."
        ),
    )
    parser.add_argument(
        "--flank-length",
        type=int,
        default=1000,
        help="Base flank length on each side of the QTL. Default: 1000.",
    )
    parser.add_argument(
        "--flank-jitter",
        type=int,
        default=300,
        help=(
            "Maximum number of bases to trim from the shared flanks independently "
            "for each genotype. Default: 300."
        ),
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=0,
        help="Seed for reproducible random generation. Default: 0.",
    )
    return parser.parse_args()


def normalize_sequence(sequence: str) -> str:
    """Remove whitespace and uppercase the sequence."""
    return re.sub(r"\s+", "", sequence).upper()


def wrap_sequence(sequence: str, line_width: int = FASTA_LINE_WIDTH) -> str:
    """Wrap a sequence to a fixed FASTA line width."""
    return "\n".join(
        sequence[i : i + line_width] for i in range(0, len(sequence), line_width)
    )


def read_fasta_records(fasta_path: Path) -> dict[str, str]:
    """Read a FASTA file into a name-to-sequence dictionary."""
    records: dict[str, list[str]] = {}
    current_name: str | None = None

    with fasta_path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current_name = line[1:].strip()
                records[current_name] = []
                continue
            if current_name is None:
                raise ValueError(f"Invalid FASTA file: {fasta_path}")
            records[current_name].append(line)

    return {
        name: normalize_sequence("".join(chunks)) for name, chunks in records.items()
    }


def write_fasta_record(record_name: str, sequence: str, output_path: Path) -> None:
    """Write a single FASTA record to a file."""
    with output_path.open("w", encoding="utf-8") as handle:
        handle.write(f">{record_name}\n")
        handle.write(f"{wrap_sequence(sequence)}\n")


def write_fasta_records(
    records: list[tuple[str, str]],
    output_path: Path,
) -> None:
    """Write multiple FASTA records to a file."""
    with output_path.open("w", encoding="utf-8") as handle:
        for record_name, sequence in records:
            handle.write(f">{record_name}\n")
            handle.write(f"{wrap_sequence(sequence)}\n")


def make_synthetic_homoeologue(sequence: str) -> str:
    """Create a highly similar same-length homoeologue by substitutions only."""
    substitutions = {"A": "C", "C": "G", "G": "T", "T": "A"}
    bases = list(sequence)

    start_index = SYNTHETIC_HOMOEOLOGUE_MUTATION_OFFSET - 1
    for index in range(
        start_index,
        len(bases),
        SYNTHETIC_HOMOEOLOGUE_MUTATION_INTERVAL,
    ):
        bases[index] = substitutions[bases[index]]

    return "".join(bases)


@define(frozen=True)
class SpacerLengthSpec:
    """Represent either a fixed spacer length or a random length range."""

    min_length: int = field(validator=instance_of(int))
    max_length: int = field(validator=instance_of(int))

    def __attrs_post_init__(self) -> None:
        """Validate spacer length bounds."""
        if self.min_length < 0 or self.max_length < 0:
            raise ValueError("Spacer lengths must be >= 0")
        if self.min_length > self.max_length:
            raise ValueError("min_length cannot be greater than max_length")

    @classmethod
    def from_string(cls, value: str) -> "SpacerLengthSpec":
        """Parse a spacer length specification."""
        cleaned = value.strip()

        if re.fullmatch(r"\d+", cleaned):
            length = int(cleaned)
            return cls(min_length=length, max_length=length)

        match = re.fullmatch(r"(\d+)-(\d+)", cleaned)
        if match:
            return cls(
                min_length=int(match.group(1)),
                max_length=int(match.group(2)),
            )

        raise ValueError("Invalid --spacer-length value. Use '300' or '250-350'.")

    def draw(self, rng: random.Random) -> int:
        """Draw one spacer length."""
        return rng.randint(self.min_length, self.max_length)


@define(frozen=True)
class Segment:
    """Represent one contiguous sequence segment."""

    segment_type: str = field(validator=instance_of(str))
    segment_name: str = field(validator=instance_of(str))
    sequence: str = field(validator=instance_of(str))

    @property
    def length(self) -> int:
        """Return segment length."""
        return len(self.sequence)


@define(frozen=True)
class SequenceAssembly:
    """Store one assembled sequence and its ordered segments."""

    sequence: str = field(validator=instance_of(str))
    segments: list[Segment] = field(validator=instance_of(list))


@define(frozen=True)
class Block:
    """Represent one named block across all genotypes."""

    name: str = field(validator=instance_of(str))
    sequences: dict[str, str] = field(validator=instance_of(dict))

    @classmethod
    def from_fasta(cls, fasta_path: Path) -> "Block":
        """Create a block from a FASTA file."""
        return cls(
            name=fasta_path.stem,
            sequences=read_fasta_records(fasta_path),
        )

    def validate_genotype_names(self, genotype_names: list[str]) -> None:
        """Check that the FASTA contains exactly the expected genotype names."""
        observed_names = sorted(self.sequences.keys())
        expected_names = sorted(genotype_names)
        if observed_names != expected_names:
            raise ValueError(
                f"Unexpected sequence names in block '{self.name}'. "
                f"Expected {expected_names}, got {observed_names}."
            )

    def get_sequence(self, genotype_name: str) -> str:
        """Return the sequence for one genotype."""
        return self.sequences[genotype_name]


@define
class BlockCollection:
    """Store all blocks loaded from the input directory."""

    blocks: dict[str, Block] = field(factory=dict, validator=instance_of(dict))

    @classmethod
    def from_directory(cls, input_dir: Path) -> "BlockCollection":
        """Load all FASTA files from a directory."""
        if not input_dir.is_dir():
            raise ValueError(f"Invalid input directory: {input_dir}")

        blocks: dict[str, Block] = {}
        for fasta_path in sorted(input_dir.iterdir()):
            if fasta_path.suffix.lower() not in VALID_FASTA_EXTENSIONS:
                continue
            block = Block.from_fasta(fasta_path)
            blocks[block.name] = block

        return cls(blocks=blocks)

    def validate_requested_blocks(self, block_order: list[str]) -> None:
        """Check that all requested blocks exist."""
        missing_blocks = [
            block_name for block_name in block_order if block_name not in self.blocks
        ]
        if missing_blocks:
            raise ValueError(f"Missing block files for: {missing_blocks}")

    def validate_genotype_names(self, genotype_names: list[str]) -> None:
        """Validate genotype names in requested blocks."""
        for block_name in BLOCK_ORDER:
            self.blocks[block_name].validate_genotype_names(genotype_names)

    def get_block(self, block_name: str) -> Block:
        """Return a block by name."""
        return self.blocks[block_name]


@define
class SequenceSpacer:
    """Generate random inter-block and edge spacer sequences."""

    rng: random.Random = field(validator=instance_of(random.Random))
    length_spec: SpacerLengthSpec = field(validator=instance_of(SpacerLengthSpec))

    def random_dna(self, length: int) -> str:
        """Generate a random DNA sequence."""
        return "".join(self.rng.choices(DNA_ALPHABET, k=length))

    def make_spacers(self, genotype_names: list[str]) -> dict[str, str]:
        """Generate one independent spacer per genotype."""
        return {
            genotype_name: self.random_dna(self.length_spec.draw(self.rng))
            for genotype_name in genotype_names
        }


@define
class SyntheticRegionBuilder:
    """Concatenate ordered blocks and spacers into one QTL per genotype."""

    block_collection: BlockCollection = field(validator=instance_of(BlockCollection))
    genotype_names: list[str] = field(validator=instance_of(list))
    spacer_generator: SequenceSpacer = field(validator=instance_of(SequenceSpacer))

    def build_qtl_assemblies(
        self, block_order: list[str]
    ) -> dict[str, SequenceAssembly]:
        """Build one QTL assembly per genotype."""
        assemblies: dict[str, list[Segment]] = {
            genotype_name: [] for genotype_name in self.genotype_names
        }

        spacer_counter = 1

        # Left edge spacer
        left_edge_spacers = self.spacer_generator.make_spacers(self.genotype_names)
        for genotype_name in self.genotype_names:
            assemblies[genotype_name].append(
                Segment(
                    segment_type="spacer",
                    segment_name=f"spacer_{spacer_counter}",
                    sequence=left_edge_spacers[genotype_name],
                )
            )
        spacer_counter += 1

        for block_index, block_name in enumerate(block_order):
            block = self.block_collection.get_block(block_name)

            # Add block
            for genotype_name in self.genotype_names:
                assemblies[genotype_name].append(
                    Segment(
                        segment_type="block",
                        segment_name=block_name,
                        sequence=block.get_sequence(genotype_name),
                    )
                )

            # Skip spacer after last block
            if block_index == len(block_order) - 1:
                continue

            # Inter-block spacer
            inter_block_spacers = self.spacer_generator.make_spacers(
                self.genotype_names
            )
            for genotype_name in self.genotype_names:
                assemblies[genotype_name].append(
                    Segment(
                        segment_type="spacer",
                        segment_name=f"spacer_{spacer_counter}",
                        sequence=inter_block_spacers[genotype_name],
                    )
                )
            spacer_counter += 1

        # Right edge spacer
        right_edge_spacers = self.spacer_generator.make_spacers(self.genotype_names)
        for genotype_name in self.genotype_names:
            assemblies[genotype_name].append(
                Segment(
                    segment_type="spacer",
                    segment_name=f"spacer_{spacer_counter}",
                    sequence=right_edge_spacers[genotype_name],
                )
            )

        return {
            genotype_name: SequenceAssembly(
                sequence="".join(segment.sequence for segment in segments),
                segments=segments,
            )
            for genotype_name, segments in assemblies.items()
        }


@define
class FlankGenerator:
    """Generate shared pseudo-chromosome flanks trimmed per genotype."""

    rng: random.Random = field(validator=instance_of(random.Random))
    flank_length: int = field(validator=instance_of(int))
    flank_jitter: int = field(validator=instance_of(int))

    def __attrs_post_init__(self) -> None:
        """Validate flank settings."""
        if self.flank_length < 0:
            raise ValueError("flank_length must be >= 0")
        if self.flank_jitter < 0:
            raise ValueError("flank_jitter must be >= 0")

    def random_dna(self, length: int) -> str:
        """Generate a random DNA sequence."""
        return "".join(self.rng.choices(DNA_ALPHABET, k=length))

    def generate_shared_flanks(self) -> tuple[str, str]:
        """Generate shared left and right flank master sequences."""
        master_length = self.flank_length + self.flank_jitter
        left_flank = self.random_dna(master_length)
        right_flank = self.random_dna(master_length)
        return left_flank, right_flank

    def trim_left_flank(self, left_flank: str) -> str:
        """Trim the shared left flank for one genotype."""
        trim_bases = self.rng.randint(0, self.flank_jitter)
        return left_flank[trim_bases:]

    def trim_right_flank(self, right_flank: str) -> str:
        """Trim the shared right flank for one genotype."""
        trim_bases = self.rng.randint(0, self.flank_jitter)
        if trim_bases == 0:
            return right_flank
        return right_flank[:-trim_bases]


@define(frozen=True)
class PseudoChromosomeRecord:
    """Store one genotype-specific QTL and pseudo-chromosome assemblies."""

    genotype_name: str = field(validator=instance_of(str))
    qtl_assembly: SequenceAssembly = field(validator=instance_of(SequenceAssembly))
    pseudochromosome_assembly: SequenceAssembly = field(
        validator=instance_of(SequenceAssembly)
    )
    qtl_start_1based: int = field(validator=instance_of(int))


@define
class DatasetBuilder:
    """Build dataset records for all genotypes."""

    flank_generator: FlankGenerator = field(validator=instance_of(FlankGenerator))

    def build_records(
        self,
        qtl_assemblies: dict[str, SequenceAssembly],
        genotype_names: list[str],
    ) -> dict[str, PseudoChromosomeRecord]:
        """Build pseudo-chromosome records for all genotypes."""
        left_flank, right_flank = self.flank_generator.generate_shared_flanks()
        records: dict[str, PseudoChromosomeRecord] = {}

        for genotype_name in genotype_names:
            trimmed_left_flank = self.flank_generator.trim_left_flank(left_flank)
            trimmed_right_flank = self.flank_generator.trim_right_flank(right_flank)
            qtl_assembly = qtl_assemblies[genotype_name]

            pseudo_segments = [
                Segment(
                    segment_type="left_flank",
                    segment_name="shared_left_flank",
                    sequence=trimmed_left_flank,
                ),
                *qtl_assembly.segments,
                Segment(
                    segment_type="right_flank",
                    segment_name="shared_right_flank",
                    sequence=trimmed_right_flank,
                ),
            ]

            pseudo_sequence = "".join(segment.sequence for segment in pseudo_segments)
            qtl_start_1based = len(trimmed_left_flank) + 1

            records[genotype_name] = PseudoChromosomeRecord(
                genotype_name=genotype_name,
                qtl_assembly=qtl_assembly,
                pseudochromosome_assembly=SequenceAssembly(
                    sequence=pseudo_sequence,
                    segments=pseudo_segments,
                ),
                qtl_start_1based=qtl_start_1based,
            )

        return records


@define
class DatasetWriter:
    """Write dataset outputs to disk."""

    output_dir: Path = field(validator=instance_of(Path))

    def fasta_record_name(self, genotype_name: str) -> str:
        """Return the target pseudochromosome FASTA record name."""
        return FASTA_HEADER_ALIASES.get(genotype_name, genotype_name)

    def homoeologue_record_name(self, genotype_name: str) -> str:
        """Return the synthetic homoeologue FASTA record name."""
        return (
            f"{self.fasta_record_name(genotype_name)}"
            f"{SYNTHETIC_HOMOEOLOGUE_SUFFIX}"
        )

    def qtl_dir(self) -> Path:
        """Return the QTL FASTA output directory."""
        return self.output_dir / "qtl_fastas"

    def pseudochromosome_dir(self) -> Path:
        """Return the pseudo-chromosome FASTA output directory."""
        return self.output_dir / "pseudochromosomes"

    def genotypes_tsv_path(self) -> Path:
        """Return the public genotypes.tsv path."""
        return self.output_dir / "genotypes.tsv"

    def chromosomes_tsv_path(self) -> Path:
        """Return the public chromosomes.tsv path."""
        return self.output_dir / "chromosomes.tsv"

    def qtl_summary_tsv_path(self) -> Path:
        """Return the QTL summary TSV path."""
        return self.output_dir / "qtl_summary.tsv"

    def pseudochromosome_summary_tsv_path(self) -> Path:
        """Return the pseudochromosome summary TSV path."""
        return self.output_dir / "pseudochromosome_summary.tsv"

    def create_directories(self) -> None:
        """Create output directories."""
        self.qtl_dir().mkdir(parents=True, exist_ok=True)
        self.pseudochromosome_dir().mkdir(parents=True, exist_ok=True)

    def qtl_fasta_path(self, genotype_name: str) -> Path:
        """Return the QTL FASTA path for one genotype."""
        return self.qtl_dir() / f"{genotype_name}.fasta"

    def pseudochromosome_fasta_path(self, genotype_name: str) -> Path:
        """Return the pseudo-chromosome FASTA path for one genotype."""
        return self.pseudochromosome_dir() / f"{genotype_name}.fasta"

    def write_fastas(self, records: dict[str, PseudoChromosomeRecord]) -> None:
        """Write regional FASTAs and two-record synthetic genome FASTAs."""
        for genotype_name, record in records.items():
            target_record_name = self.fasta_record_name(genotype_name)
            homoeologue_record_name = self.homoeologue_record_name(genotype_name)

            write_fasta_record(
                record_name=target_record_name,
                sequence=record.qtl_assembly.sequence,
                output_path=self.qtl_fasta_path(genotype_name),
            )

            target_genome_sequence = record.pseudochromosome_assembly.sequence
            homoeologue_sequence = make_synthetic_homoeologue(
                target_genome_sequence
            )
            write_fasta_records(
                records=[
                    (target_record_name, target_genome_sequence),
                    (homoeologue_record_name, homoeologue_sequence),
                ],
                output_path=self.pseudochromosome_fasta_path(genotype_name),
            )

    def write_genotypes_tsv(self, records: dict[str, PseudoChromosomeRecord]) -> None:
        """Write the headered genotypes.tsv file."""
        with self.genotypes_tsv_path().open("w", encoding="utf-8") as handle:
            handle.write("genotype\tgroup\tregion_fasta\tgenome_fasta\tsource_seq\tregion_start\n")
            for genotype_name in GENOTYPE_NAMES:
                record = records[genotype_name]
                qtl_fasta_path = self.qtl_fasta_path(genotype_name)
                genome_fasta_path = self.pseudochromosome_fasta_path(genotype_name)
                group = {"Karur": "A", "Soldur": "B"}.get(genotype_name, "")
                handle.write(
                    f"{genotype_name}\t{group}\t{qtl_fasta_path}\t{genome_fasta_path}\t"
                    f"{self.fasta_record_name(genotype_name)}\t{record.qtl_start_1based}\n"
                )

    def write_chromosomes_tsv(
        self,
        records: dict[str, PseudoChromosomeRecord],
    ) -> None:
        """Write chromosome/subgenome metadata for generated genome FASTAs."""
        with self.chromosomes_tsv_path().open("w", encoding="utf-8") as handle:
            handle.write(
                "genotype\tseq_id\thomoeologous_group\tsubgenome\n"
            )
            for genotype_name in GENOTYPE_NAMES:
                # Access the record so missing genotypes fail consistently.
                records[genotype_name]

                handle.write(
                    f"{genotype_name}\t{self.fasta_record_name(genotype_name)}\t"
                    f"{SYNTHETIC_HOMOEOLOGOUS_GROUP}\t"
                    f"{SYNTHETIC_TARGET_SUBGENOME}\n"
                )
                handle.write(
                    f"{genotype_name}\t{self.homoeologue_record_name(genotype_name)}\t"
                    f"{SYNTHETIC_HOMOEOLOGOUS_GROUP}\t"
                    f"{SYNTHETIC_HOMOEOLOGUE_SUBGENOME}\n"
                )

    def write_summary_tsvs(self, records: dict[str, PseudoChromosomeRecord]) -> None:
        """Write separate QTL and pseudochromosome summary TSV files."""
        with self.qtl_summary_tsv_path().open("w", encoding="utf-8") as qtl_handle:
            qtl_handle.write("genotype\tsegment\tstart_1based\tend_1based\tlength\n")
            for genotype_name in GENOTYPE_NAMES:
                record = records[genotype_name]
                self.write_assembly_summary(
                    handle=qtl_handle,
                    genotype_name=genotype_name,
                    assembly=record.qtl_assembly,
                )

        with self.pseudochromosome_summary_tsv_path().open(
            "w",
            encoding="utf-8",
        ) as pseudo_handle:
            pseudo_handle.write("genotype\tsegment\tstart_1based\tend_1based\tlength\n")
            for genotype_name in GENOTYPE_NAMES:
                record = records[genotype_name]
                self.write_assembly_summary(
                    handle=pseudo_handle,
                    genotype_name=genotype_name,
                    assembly=record.pseudochromosome_assembly,
                )

    @staticmethod
    def write_assembly_summary(
        handle,
        genotype_name: str,
        assembly: SequenceAssembly,
    ) -> None:
        """Write one assembly summary to an open TSV handle."""
        start_1based = 1

        for segment in assembly.segments:
            end_1based = start_1based + segment.length - 1
            handle.write(
                f"{genotype_name}\t{segment.segment_name}\t"
                f"{start_1based}\t{end_1based}\t{segment.length}\n"
            )
            start_1based = end_1based + 1


def main() -> None:
    """Run the program."""
    args = parse_args()

    rng = random.Random(args.seed)
    spacer_length_spec = SpacerLengthSpec.from_string(args.spacer_length)

    block_collection = BlockCollection.from_directory(Path(args.input_dir))
    block_collection.validate_requested_blocks(BLOCK_ORDER)
    block_collection.validate_genotype_names(GENOTYPE_NAMES)

    spacer_generator = SequenceSpacer(
        rng=rng,
        length_spec=spacer_length_spec,
    )
    synthetic_region_builder = SyntheticRegionBuilder(
        block_collection=block_collection,
        genotype_names=GENOTYPE_NAMES,
        spacer_generator=spacer_generator,
    )
    qtl_assemblies = synthetic_region_builder.build_qtl_assemblies(BLOCK_ORDER)

    flank_generator = FlankGenerator(
        rng=rng,
        flank_length=args.flank_length,
        flank_jitter=args.flank_jitter,
    )
    dataset_builder = DatasetBuilder(flank_generator=flank_generator)
    records = dataset_builder.build_records(qtl_assemblies, GENOTYPE_NAMES)

    dataset_writer = DatasetWriter(output_dir=Path(args.output_dir))
    dataset_writer.create_directories()
    dataset_writer.write_fastas(records)
    dataset_writer.write_genotypes_tsv(records)
    dataset_writer.write_chromosomes_tsv(records)
    dataset_writer.write_summary_tsvs(records)


if __name__ == "__main__":
    main()
