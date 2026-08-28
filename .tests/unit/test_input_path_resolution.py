from __future__ import annotations

import copy
import os
import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "workflow" / "scripts"))

from run_info import start_run_record
from validate_inputs import validate_inputs


def create_analysis(analysis_dir: Path) -> dict:
    """Create a minimal self-contained KASP analysis directory."""
    (analysis_dir / "regions").mkdir(parents=True)
    (analysis_dir / "genomes").mkdir()
    (analysis_dir / "annotations").mkdir()

    (analysis_dir / "regions" / "A.fasta").write_text(
        ">region\nACGTACGT\n",
        encoding="utf-8",
    )
    (analysis_dir / "genomes" / "A.fasta").write_text(
        ">chr1\nACGTACGT\n",
        encoding="utf-8",
    )
    (analysis_dir / "annotations" / "A.gff").write_text(
        "chr1\tsource\tgene\t1\t8\t.\t+\t.\tID=gene1\n",
        encoding="utf-8",
    )
    (analysis_dir / "repeats.fa").write_text(
        ">repeat\nACGT\n",
        encoding="utf-8",
    )
    (analysis_dir / "genotypes.tsv").write_text(
        "genotype\tgroup\tregion_fasta\tgenome_fasta\t"
        "source_seq\tregion_start\n"
        "A\t\tregions/A.fasta\tgenomes/A.fasta\tchr1\t1\n",
        encoding="utf-8",
    )
    (analysis_dir / "annotations.tsv").write_text(
        "genotype\ttrack\tgff\nA\tgenes\tannotations/A.gff\n",
        encoding="utf-8",
    )
    (analysis_dir / "chromosomes.tsv").write_text(
        "genotype\tseq_id\thomoeologous_group\tsubgenome\n"
        "A\tchr1\t1\tA\n",
        encoding="utf-8",
    )

    config = {
        "project": {"name": "Path resolution test"},
        "inputs": {
            "genotypes": "genotypes.tsv",
            "annotations": "annotations.tsv",
            "chromosomes": "chromosomes.tsv",
        },
        "snps": {
            "repeat_masking": {"library": "repeats.fa"},
        },
        "kasp": {"polymarker_genomes": 1},
    }
    (analysis_dir / "config.yaml").write_text(
        "project:\n"
        "  name: Path resolution test\n"
        "inputs:\n"
        "  genotypes: genotypes.tsv\n"
        "  annotations: annotations.tsv\n"
        "  chromosomes: chromosomes.tsv\n"
        "snps:\n"
        "  repeat_masking:\n"
        "    library: repeats.fa\n"
        "kasp:\n"
        "  polymarker_genomes: 1\n",
        encoding="utf-8",
    )
    return config


def validate_kasp_config(config: dict, working_directory: Path) -> None:
    """Validate a complete KASP config using the supplied workdir."""
    validate_inputs(
        config=copy.deepcopy(config),
        mode="kasp",
        schema_dir=REPO_ROOT / "workflow" / "schemas",
        working_directory=working_directory,
    )


def test_validation_uses_snakemake_working_directory(
    tmp_path: Path,
    monkeypatch,
) -> None:
    """Relative config and table paths are rooted at --directory."""
    analysis_dir = tmp_path / "analysis"
    caller_dir = tmp_path / "caller"
    caller_dir.mkdir()
    config = create_analysis(analysis_dir)

    monkeypatch.chdir(caller_dir)

    validate_kasp_config(config, analysis_dir)


def test_validation_accepts_absolute_config_paths(
    tmp_path: Path,
    monkeypatch,
) -> None:
    """Absolute config paths remain independent of the working directory."""
    analysis_dir = tmp_path / "analysis"
    caller_dir = tmp_path / "caller"
    caller_dir.mkdir()
    config = create_analysis(analysis_dir)
    config["inputs"] = {
        key: str(analysis_dir / f"{key}.tsv")
        for key in ("genotypes", "annotations", "chromosomes")
    }
    config["snps"]["repeat_masking"]["library"] = str(
        analysis_dir / "repeats.fa"
    )

    monkeypatch.chdir(caller_dir)

    validate_kasp_config(config, analysis_dir)


def test_validation_without_directory_uses_current_working_directory(
    tmp_path: Path,
    monkeypatch,
) -> None:
    """Relative paths keep their current-directory behavior without -d."""
    analysis_dir = tmp_path / "analysis"
    config = create_analysis(analysis_dir)

    monkeypatch.chdir(analysis_dir)

    validate_kasp_config(config, Path.cwd())


def test_provenance_snapshots_inputs_from_snakemake_working_directory(
    tmp_path: Path,
    monkeypatch,
) -> None:
    """Snapshots use relative input paths from the selected workdir."""
    analysis_dir = tmp_path / "analysis"
    caller_dir = tmp_path / "caller"
    caller_dir.mkdir()
    config = create_analysis(analysis_dir)

    monkeypatch.chdir(caller_dir)

    record = start_run_record(
        config=config,
        mode="kasp",
        config_path=analysis_dir / "config.yaml",
        repo_root=tmp_path,
        argv=["kasp", "-c", "config.yaml"],
        snakemake_args=["--directory", str(analysis_dir)],
        working_directory=analysis_dir,
    )

    assert record.history_path.parent == analysis_dir / "results" / "run_history"
    snapshots = record.data["inputs"]["snapshots"]
    for key in ("config_file", "genotypes", "annotations", "chromosomes"):
        snapshot_path = Path(snapshots[key])
        assert snapshot_path.is_file()
        assert snapshot_path.parent == record.history_path.with_suffix("")


def test_cli_dry_run_uses_directory_for_relative_inputs(
    tmp_path: Path,
) -> None:
    """CLI validation and Snakemake agree on relative paths with -d."""
    analysis_dir = tmp_path / "analysis"
    caller_dir = tmp_path / "caller"
    caller_dir.mkdir()
    create_analysis(analysis_dir)

    result = subprocess.run(
        [
            sys.executable,
            str(REPO_ROOT / "workflow" / "scripts" / "kaspberry.py"),
            "kasp",
            "-c",
            str(analysis_dir / "config.yaml"),
            "--directory",
            str(analysis_dir),
            "-n",
        ],
        cwd=caller_dir,
        env=os.environ.copy(),
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stderr
    assert "results/kasp/kasp_snp_summary.tsv" in result.stdout
