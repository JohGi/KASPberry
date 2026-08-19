#!/usr/bin/env python3
# Author: Johanna Girodolle

from __future__ import annotations

from pathlib import Path

from attrs import define, field

@define(frozen=True)
class SampleRecord:
    """Store one sample definition from the samples TSV."""

    fasta_path: Path
    sample: str
    region_start_in_source_seq: int = 1


@define(frozen=True)
class BlockFeature:
    """Store one collinear block for one sample."""

    sample: str
    block_id: str
    block_start_in_region: int
    block_end_in_region: int
    block_start_in_source_seq: int
    block_end_in_source_seq: int

    @property
    def feature_id(self) -> str:
        """Return a shared feature identifier."""
        return f"block::{self.block_id}"


@define(frozen=True)
class SnpFeature:
    """Store one SNP for one sample."""

    sample: str
    block_id: str
    aln_pos: int
    nt: str
    pos_in_block: int
    pos_in_region: int
    pos_in_source_seq: int

    @property
    def feature_id(self) -> str:
        """Return a shared feature identifier."""
        return f"snp::{self.block_id}::{self.aln_pos}"


@define(frozen=True)
class SampleData:
    """Store all display data for one sample."""

    sample: str
    region_length: int
    region_start_in_source_seq: int = 1
    blocks: list[BlockFeature] = field(factory=list)
    snps: list[SnpFeature] = field(factory=list)


@define(frozen=True)
class DistanceMatrix:
    """Represent a square distance matrix."""

    labels: list[str]
    values: list[list[float]]
    source: str
    title: str
    unit: str

    def to_dict(self) -> dict[str, object]:
        """Convert the distance matrix to a JSON-compatible dictionary."""
        return {
            "labels": self.labels,
            "values": self.values,
            "source": self.source,
            "title": self.title,
            "unit": self.unit,
        }


@define(frozen=True)
class BlockAlignment:
    """Store one multiple-sequence alignment for a collinear block."""

    block_id: str
    sequences_by_sample: dict[str, str]

    def __attrs_post_init__(self) -> None:
        """Validate that all aligned sequences have the same length."""
        lengths = {len(sequence) for sequence in self.sequences_by_sample.values()}

        if len(lengths) > 1:
            raise ValueError(
                f"Alignment for block {self.block_id} contains sequences "
                "with inconsistent lengths."
            )

    @property
    def length(self) -> int:
        """Return the alignment length."""
        if not self.sequences_by_sample:
            return 0

        return len(next(iter(self.sequences_by_sample.values())))

    def to_payload(self) -> dict[str, str]:
        """Return the alignment payload indexed by sample name."""
        return self.sequences_by_sample


@define(frozen=True)
class GffGeneFeature:
    """Store one GFF gene feature projected into the displayed region."""

    sample: str
    track_name: str
    gene_id: str
    source_seq_id: str
    start_in_source_seq: int
    end_in_source_seq: int
    start_in_region: int
    end_in_region: int
    strand: str | None = None

@define(frozen=True)
class GffTrack:
    """Store projected GFF gene features for one sample track."""

    sample: str
    track_name: str
    features: list[GffGeneFeature]


@define(frozen=True)
class DotplotRecord:
    """Store one dotplot entry from the dotplot manifest JSON."""

    pair_id: str
    x_sample: str
    y_sample: str
    svg_rel_path: str
