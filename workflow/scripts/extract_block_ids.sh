#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Extract unique block IDs from a GFF file.

Usage:
  extract_block_ids.sh <input.gff> <output.list>
EOF
}

if [[ $# -ne 2 ]]; then
    usage >&2
    exit 1
fi

input_gff="$1"
output_list="$2"

if [[ ! -s "$input_gff" ]]; then
    echo "Error: input GFF '$input_gff' not found or empty." >&2
    exit 1
fi

awk -F'\t' '
    $0 !~ /^#/ {
        n = split($9, attrs, ";")
        for (i = 1; i <= n; i++) {
            if (attrs[i] ~ /^ID=/) {
                sub(/^ID=/, "", attrs[i])
                print attrs[i]
            }
        }
    }
' "$input_gff" | sort -u > "$output_list"
