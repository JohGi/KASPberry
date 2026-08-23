#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Run blastn2dotplots for one pair.

Usage:
  run_pairwise_blastn2dotplots.sh \
    --blastn-tsv pair.tsv \
    --db-name Target \
    --query-name Query \
    --out-prefix out/prefix \
EOF
}

blastn_tsv=""
db_name=""
query_name=""
out_prefix=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --blastn-tsv)
            blastn_tsv="${2:-}"
            shift 2
            ;;
        --db-name)
            db_name="${2:-}"
            shift 2
            ;;
        --query-name)
            query_name="${2:-}"
            shift 2
            ;;
        --out-prefix)
            out_prefix="${2:-}"
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

if [[ -z "$blastn_tsv" || -z "$db_name" || -z "$query_name" || -z "$out_prefix" ]]; then
    echo "Error: missing required arguments." >&2
    usage >&2
    exit 1
fi

if [[ ! -s "$blastn_tsv" ]]; then
    echo "Error: BLAST-like TSV '$blastn_tsv' not found or empty." >&2
    exit 1
fi

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"

mkdir -p "$(dirname "$out_prefix")"

tmp_dir="$(mktemp -d)"
trap 'rm -rf "$tmp_dir"' EXIT

db_file="$tmp_dir/db.txt"
query_file="$tmp_dir/query.txt"

printf "%s\t%s\n" "$db_name" "$db_name" > "$db_file"
printf "%s\t%s\n" "$query_name" "$query_name" > "$query_file"

python "${SCRIPT_DIR}/blastn2dotplots/blastn2dotplots" \
    -i1 "$db_file" \
    -i2 "$query_file" \
    --blastn "$blastn_tsv" \
    --out "$out_prefix"
