# KASPberry configuration

KASPberry is configured with a YAML file plus one required genotype table and, depending on the analysis, optional annotation and chromosome tables.

The authoritative list of configuration keys, types, defaults, and validation constraints is defined in `workflow/schemas/config.schema.yaml`. Input tables are validated against the corresponding schemas in `workflow/schemas/`.

## Running KASPberry

For normal analyses, use the KASPberry command-line interface and provide the configuration explicitly:

```bash
kaspberry snps     -c config/config.snps.yaml     --directory /path/to/analysis     --cores all
```

or:

```bash
kaspberry kasp     -c config/config.kasp.yaml     --directory /path/to/analysis     --cores all
```

KASPberry forwards additional command-line arguments to Snakemake. Results are written to `results/` inside the Snakemake working directory, so the commands above write to `/path/to/analysis/results/`.

Paths in the YAML file and TSV input tables are interpreted relative to the Snakemake working directory selected with `--directory`, not relative to the location of the YAML file. Absolute paths can also be used.

Two example configuration files are provided and can be copied and adapted as starting points:

- `config/config.snps.yaml` for SNP discovery and diagnostic SNP filtering;
- `config/config.kasp.yaml` for KASP assay design and in silico validation.
  

## Configuration file

The two modes share most settings, while the KASP workflow requires additional inputs and assay-design parameters. The complete configuration schema is defined in `workflow/schemas/config.schema.yaml`.  
A typical configuration has the following structure:

```yaml
project:
  name: "My analysis"

inputs:
  genotypes: genotypes.tsv
  annotations: null
  chromosomes: null

snps:
  min_block_length: 301
  repeat_masking:
    library: null

kasp:
  polymarker_genomes: 2
  mfeprimer:
    specificity:
      min_tm: 50
    dimer:
      max_dg: -3.5

viewer:
  dotplot_reference: null

advanced:
  sibeliaz:
    extra_options: []
  batching:
    blocks_per_chunk: 2
  alignment:
    mafft_options: []
  mfeprimer:
    kmer_size: 9
    specificity_extra_options: []
    dimer_extra_options: []
    hairpin_extra_options: []
```

The `kasp` section is only required for `kaspberry kasp`.

## `genotypes.tsv`

`inputs.genotypes` is required in both modes. The table must be tab-separated and contain the following columns. Columns whose values are optional must still be present in the header; unused cells can be left empty.

| Column | Required value | Description |
| --- | --- | --- |
| `genotype` | yes | Unique genotype/accession identifier. |
| `group` | for diagnostic grouping | Diagnostic group assigned to the genotype. The table must contain either no group values or exactly two distinct groups. |
| `region_fasta` | yes | FASTA containing the homologous regional sequence analysed for this genotype. It must contain exactly one non-empty FASTA record. |
| `genome_fasta` | KASP QC only | Whole-genome FASTA used for PolyMarker design and MFEprimer validation. |
| `source_seq` | when source coordinates are needed | Identifier of the chromosome/pseudomolecule containing the regional sequence. |
| `region_start` | together with `source_seq` | 1-based start coordinate of `region_fasta` on `source_seq`. |

Example:

```text
genotype	group	region_fasta	genome_fasta	source_seq	region_start
Elite1	short	regions/Elite1.fasta	genomes/Elite1.fasta	2A	15000001
Elite2	short	regions/Elite2.fasta	genomes/Elite2.fasta	2A	14987022
Wild1	tall	regions/Wild1.fasta	genomes/Wild1.fasta	2A	15210411
Wild2	tall	regions/Wild2.fasta			
```

For diagnostic SNP discovery, assign the accessions to exactly two groups. All genotypes within one group must carry the same allele and the two groups must carry different alleles for a SNP to be retained as diagnostic.

If `genome_fasta` is provided, both `source_seq` and `region_start` are required. KASPberry checks that `source_seq` exists in the genome FASTA, that the regional interval fits within that sequence, and that the beginning of `region_fasta` matches the genome sequence at the declared coordinate.

At least one genotype must provide `genome_fasta`, `source_seq`, and `region_start` in `kasp` mode.

## `annotations.tsv`

`inputs.annotations` is optional in both modes. When provided, it is a tab-separated table with the following columns:

| Column | Description |
| --- | --- |
| `genotype` | Genotype/accession to which the track belongs. |
| `track` | Display name used for the annotation track. |
| `gff` | Path to the GFF/GFF3 file. |

Example:

```text
genotype	track	gff
Elite1	Genes	annotations/Elite1.gff3
Wild1	Genes	annotations/Wild1.gff3
```

A genotype with an annotation track must also provide `source_seq` and `region_start` in `genotypes.tsv`. The declared `source_seq` must occur as a sequence identifier in the corresponding GFF/GFF3 file.

## `chromosomes.tsv`

`inputs.chromosomes` describes chromosome or pseudomolecule records whose homoeologous group and subgenome assignment is known.

It is required in `kasp` mode when `kasp.polymarker_genomes > 1` and optional when `kasp.polymarker_genomes == 1`. When chromosome or pseudomolecule assignments are known, providing this table is recommended even when `kasp.polymarker_genomes == 1`.

All records from each `genome_fasta` are retained for PolyMarker design and MFEprimer validation. When `chromosomes.tsv` is provided, records listed in the table are assigned PolyMarker-compatible identifiers preserving their declared homoeologous-group and subgenome structure. Records not listed in the table are retained but treated as unassigned sequences belonging to a single artificial PolyMarker chromosome group.

When `kasp.polymarker_genomes == 1` and no `chromosomes.tsv` is provided, each FASTA record is treated as an independent PolyMarker chromosome group.

The table must contain:

| Column               | Description                                                      |
| -------------------- | ---------------------------------------------------------------- |
| `genotype`           | Genotype/accession declared in `genotypes.tsv`.                  |
| `seq_id`             | Exact FASTA record identifier in that genotype's `genome_fasta`. |
| `homoeologous_group` | User-defined label for the homoeologous chromosome group.        |
| `subgenome`          | User-defined subgenome label.                                    |

Example for a tetraploid genome:

```text
genotype	seq_id	homoeologous_group	subgenome
Elite1	2A	2	A
Elite1	2B	2	B
Wild1	2A	2	A
Wild1	2B	2	B
```

Users are encouraged to list all known chromosomes or pseudomolecules of the assembly. Smaller contigs, scaffolds, or other unassigned sequences do not need to be listed and remain included automatically as unassigned sequences.

For every genotype used for KASP whole-genome QC, when `chromosomes.tsv` is provided, the row corresponding to `source_seq` must be present. Within its homoeologous group, the number of distinct `subgenome` labels must equal `kasp.polymarker_genomes`.

PolyMarker-compatible chromosome-group identifiers are limited to 68 distinct values. When `chromosomes.tsv` is provided and unassigned FASTA records are present, these records collectively require one additional artificial group. Therefore, at most 67 annotated homoeologous groups can coexist with unassigned sequences, or 68 annotated groups when every FASTA record is assigned. When `kasp.polymarker_genomes == 1` and no `chromosomes.tsv` is provided, the genome FASTA can therefore contain at most 68 records.

## SNP discovery settings

### `snps.min_block_length`

Minimum length, in base pairs, of conserved/collinear blocks retained for SNP discovery. Default: 301; minimum allowed value: 301. This threshold is consistent with KASPberry's fixed requirement of 150 bp of polymorphism-free aligned sequence on each side of a candidate SNP.

### `snps.repeat_masking.library`

Simple repeats and low-complexity regions are always masked before SNP discovery.

Set `library` to a FASTA repeat library to additionally mask user-defined repeat sequences, or to `null` to use only the default simple-repeat/low-complexity masking.

## KASP assay settings

### `kasp.polymarker_genomes`

Number of subgenomes expected by PolyMarker for the target homoeologous group.
For example, use `2` for a tetraploid genome when the target group contains
A- and B-subgenome chromosomes.

Allowed values range from `1` to `62`. When the value is greater than `1`,
`inputs.chromosomes` is required so that KASPberry can identify the expected
homoeologous chromosome structure.

### `kasp.mfeprimer.specificity.min_tm`

Minimum primer-binding melting temperature used during MFEprimer specificity screening. Default: `50`.

### `kasp.mfeprimer.dimer.max_dg`

Delta-G threshold used by MFEprimer for reporting potentially problematic primer dimers. Default: `-3.5`.

## Viewer settings

### `viewer.dotplot_reference`

Optional genotype to use as the common reference for pairwise dotplots. Set to `null` to use the workflow default behavior.

## Advanced settings

The `advanced` section exposes a limited number of tool-level and batching controls:

- `advanced.sibeliaz.extra_options`: additional SibeliaZ options.
- `advanced.batching.blocks_per_chunk`: number of conserved blocks processed per workflow chunk.
- `advanced.alignment.mafft_options`: additional MAFFT options.
- `advanced.mfeprimer.kmer_size`: 3' seed length used when building MFEprimer genome indexes; allowed range `1-15`, default `9`.
- `advanced.mfeprimer.specificity_extra_options`: additional `mfeprimer spec` options.
- `advanced.mfeprimer.dimer_extra_options`: additional `mfeprimer dimer` options.
- `advanced.mfeprimer.hairpin_extra_options`: additional `mfeprimer hairpin` options.

Options that would override command-line arguments managed directly by KASPberry are rejected during input validation.

## Snakemake Workflow Catalog

The workflow's default Snakemake target is the `snps` mode. Direct Snakefile invocation is supported for workflow inspection and Snakemake Workflow Catalog integration, but the KASPberry CLI is the supported interface for selecting between `snps` and `kasp` in normal analyses.

The `.test/` directory contains the self-contained smoke-test dataset used by the Catalog to construct the workflow rule graph. It is not intended as a user configuration.
