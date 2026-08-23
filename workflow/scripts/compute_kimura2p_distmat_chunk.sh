#!/usr/bin/env bash
# Author: Johanna Girodolle

set -euo pipefail

usage() {
    cat <<'EOF'
Run EMBOSS distmat with the Kimura 2-parameter nucleotide model for one chunk.

Usage:
  compute_kimura2p_distmat_chunk.sh \
    --chunk-list chunk_0000.list \
    --aln-dir results/15_unmasked_alignments \
    --outdir results/16_kimura2p_distances/matrices

Description:
  Compute one EMBOSS distmat matrix per aligned unmasked block FASTA.
  The distance method is fixed to nucmethod 2, corresponding to Kimura 2-parameter distance.

  FASTA headers are cleaned before running EMBOSS:
    >Belalur:3359-3893 becomes >Belalur
EOF
}

chunk_list=""
aln_dir=""
outdir=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --chunk-list)
            chunk_list="${2:-}"
            shift 2
            ;;
        --aln-dir)
            aln_dir="${2:-}"
            shift 2
            ;;
        --outdir)
            outdir="${2:-}"
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

if [[ -z "$chunk_list" || -z "$aln_dir" || -z "$outdir" ]]; then
    echo "Error: missing required arguments." >&2
    usage >&2
    exit 1
fi

if [[ ! -f "$chunk_list" || ! -s "$chunk_list" ]]; then
    echo "Error: chunk list '$chunk_list' not found, not a file, or empty." >&2
    exit 1
fi

if [[ ! -d "$aln_dir" ]]; then
    echo "Error: alignment directory '$aln_dir' not found." >&2
    exit 1
fi

tmp_dir="$(mktemp -d)"
trap 'rm -rf "$tmp_dir"' EXIT

mkdir -p "$outdir"

clean_fasta_headers() {
    local input_fasta="$1"
    local output_fasta="$2"

    awk '
        /^>/ {
            header = substr($0, 2)
            sub(/:[^:]*$/, "", header)
            print ">" header
            next
        }
        {
            print
        }
    ' "$input_fasta" > "$output_fasta"
}

while read -r block_id; do
    [[ -z "$block_id" ]] && continue

    aln_path="${aln_dir}/${block_id}.aln.fasta"
    cleaned_aln_path="${tmp_dir}/${block_id}.cleaned.aln.fasta"
    output_path="${outdir}/${block_id}.kimura2p.distmat"

    if [[ ! -s "$aln_path" ]]; then
        echo "Error: alignment '$aln_path' not found or empty." >&2
        exit 1
    fi

    echo "[INFO] Cleaning FASTA headers for block '${block_id}'." >&2
    clean_fasta_headers "$aln_path" "$cleaned_aln_path"

    echo "[INFO] Computing Kimura 2-parameter distance matrix for block '${block_id}'." >&2
    distmat \
        -sequence "$cleaned_aln_path" \
        -outfile "$output_path" \
        -nucmethod 2
done < "$chunk_list"
