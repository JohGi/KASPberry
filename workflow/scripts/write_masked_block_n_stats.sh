#!/usr/bin/env bash
# Author: Johanna Girodolle

set -euo pipefail

usage() {
    cat <<EOF
Usage:
  write_masked_block_n_stats.sh \
    --masked-dir <dir> \
    --unmasked-dir <dir> \
    --output <path>

Description:
  Aggregate paired N-content statistics for all masked and unmasked block
  alignment FASTA files in the given directories.
EOF
}

parse_args() {
    MASKED_DIR=""
    UNMASKED_DIR=""
    OUTPUT=""

    while [[ $# -gt 0 ]]; do
        case "$1" in
            --masked-dir)
                MASKED_DIR="$2"
                shift 2
                ;;
            --unmasked-dir)
                UNMASKED_DIR="$2"
                shift 2
                ;;
            --output)
                OUTPUT="$2"
                shift 2
                ;;
            --help|-h)
                usage
                exit 0
                ;;
            *)
                echo "[ERROR] Unknown argument: $1" >&2
                usage >&2
                exit 1
                ;;
        esac
    done

    if [[ -z "${MASKED_DIR}" || -z "${UNMASKED_DIR}" || -z "${OUTPUT}" ]]; then
        echo "[ERROR] Missing required arguments." >&2
        usage >&2
        exit 1
    fi
}

get_python_script_path() {
    local script_dir
    script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    PYTHON_SCRIPT="${script_dir}/compute_n_stats.py"
}

write_header() {
    local output_path="$1"
    printf "block_id\tsample\tmasked_length_bp\tunmasked_length_bp\tunmasked_n_count\tunmasked_n_pct\tmasked_n_count\tmasked_n_pct\trepeat_masked_n_count\trepeat_masked_n_pct\n" \
        > "$output_path"
}

append_stats_for_all_fastas() {
    local masked_dir="$1"
    local unmasked_dir="$2"
    local python_script="$3"
    local output_path="$4"

    local block_id
    local masked_fasta
    local unmasked_fasta

    shopt -s nullglob
    for masked_fasta in "${masked_dir}"/*.aln.fasta; do
        block_id="$(basename "${masked_fasta}" .aln.fasta)"
        unmasked_fasta="${unmasked_dir}/${block_id}.aln.fasta"

        if [[ ! -f "${unmasked_fasta}" ]]; then
            echo "[ERROR] Unmasked alignment not found for block ${block_id}: ${unmasked_fasta}" >&2
            exit 1
        fi

        python3 "${python_script}" \
            --masked-alignment "${masked_fasta}" \
            --unmasked-alignment "${unmasked_fasta}" \
            >> "$output_path"
    done
    shopt -u nullglob
}

main() {
    parse_args "$@"
    get_python_script_path

    local tmp_output
    tmp_output="${OUTPUT}.tmp"

    write_header "$tmp_output"
    append_stats_for_all_fastas "$MASKED_DIR" "$UNMASKED_DIR" "$PYTHON_SCRIPT" "$tmp_output"

    mv "$tmp_output" "$OUTPUT"
}

main "$@"
