#!/usr/bin/env bash
# Author: Johanna Girodolle

set -euo pipefail

usage() {
    cat <<'EOF'
Run MAFFT sequentially for all block FASTA files listed in one chunk.

Usage:
  align_block_chunk.sh \
    --chunk-list chunk_0000.list \
    --fasta-dir results/05_masked_block_fastas \
    --outdir results/06_alignments \
    --threads 1 \
    --fasta-suffix .fasta.masked \
    --mafft-extra-options "--auto"
EOF
}

chunk_list=""
fasta_dir=""
outdir=""
threads=""
fasta_suffix=""
mafft_extra_options=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --chunk-list)
            chunk_list="${2:-}"
            shift 2
            ;;
        --fasta-dir)
            fasta_dir="${2:-}"
            shift 2
            ;;
        --outdir)
            outdir="${2:-}"
            shift 2
            ;;
        --threads)
            threads="${2:-}"
            shift 2
            ;;
        --fasta-suffix)
            fasta_suffix="${2:-}"
            shift 2
            ;;
        --mafft-extra-options)
            mafft_extra_options="${2:-}"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Error: unknown argument '$1'." >&2
            usage >&2
            exit 1
            ;;
    esac
done

if [[ -z "$chunk_list" || -z "$fasta_dir" || -z "$outdir" || -z "$threads" || -z "$fasta_suffix" ]]; then
    echo "Error: missing required arguments." >&2
    usage >&2
    exit 1
fi

if [[ ! -s "$chunk_list" ]]; then
    echo "Error: chunk list '$chunk_list' not found or empty." >&2
    exit 1
fi

if [[ ! -d "$fasta_dir" ]]; then
    echo "Error: FASTA directory '$fasta_dir' not found." >&2
    exit 1
fi

mkdir -p "$outdir"

while read -r block_id; do
    [[ -z "$block_id" ]] && continue

    fasta_path="${fasta_dir}/${block_id}${fasta_suffix}"

    if [[ ! -s "$fasta_path" ]]; then
        echo "Error: FASTA '$fasta_path' not found or empty." >&2
        exit 1
    fi

    output_path="${outdir}/${block_id}.aln.fasta"

    echo "[INFO] Aligning block '${block_id}' from '${fasta_path}'." >&2
    mafft --thread "$threads" $mafft_extra_options "$fasta_path" > "$output_path"
done < "$chunk_list"
