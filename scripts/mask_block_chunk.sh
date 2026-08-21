#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Run repeat masking sequentially for all block FASTA files listed in one chunk.

Usage:
  mask_block_chunk.sh \
    --chunk-list chunk_0000.list \
    --fasta-dir results/04_block_fastas/per_block \
    [--te-lib repeats.fasta] \
    --outdir results/05_masked_block_fastas \
    --threads 1
EOF
}

chunk_list=""
fasta_dir=""
te_lib=""
outdir=""
threads=""

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
        --te-lib)
            te_lib="${2:-}"
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

if [[ -z "$chunk_list" || -z "$fasta_dir" || -z "$outdir" || -z "$threads" ]]; then
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

if [[ -n "$te_lib" && ! -s "$te_lib" ]]; then
    echo "Error: TE library '$te_lib' not found or empty." >&2
    exit 1
fi

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mask_script="${script_dir}/mask_repeats.sh"

if [[ ! -f "$mask_script" ]]; then
    echo "Error: mask script '$mask_script' not found." >&2
    exit 1
fi

mkdir -p "$outdir"

while read -r block_id; do
    [[ -z "$block_id" ]] && continue

    fasta_path="${fasta_dir}/${block_id}.fasta"
    output_path="${outdir}/${block_id}.fasta.masked"

    echo "[INFO] Masking block '${block_id}'." >&2

    mask_args=(
        --fasta "$fasta_path"
        --outdir "$outdir"
        --threads "$threads"
        --output "$output_path"
    )

    if [[ -n "$te_lib" ]]; then
        mask_args+=( --te-lib "$te_lib" )
    fi

    bash "$mask_script" "${mask_args[@]}"
done < "$chunk_list"
