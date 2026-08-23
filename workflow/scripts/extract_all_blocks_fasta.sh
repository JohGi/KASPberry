#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Extract all block sequences from a multi-FASTA and a filtered GFF file.

Usage:
  extract_all_blocks_fasta.sh \
    --fasta all_genomes.fasta \
    --gff filtered_blocks.gff \
    --output all_blocks.raw.fasta

Notes:
  - The FASTA sequence names must match column 1 of the GFF.
  - GFF coordinates are converted to BED coordinates.
  - Strand is ignored on purpose.
  - Output FASTA headers are formatted as:
      >BLOCKID__SEQNAME:START-END
EOF
}

fasta=""
gff=""
output=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --fasta)
            fasta="${2:-}"
            shift 2
            ;;
        --gff)
            gff="${2:-}"
            shift 2
            ;;
        --output)
            output="${2:-}"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "ERROR: unknown argument '$1'." >&2
            usage >&2
            exit 1
            ;;
    esac
done

if [[ -z "$fasta" || -z "$gff" || -z "$output" ]]; then
    echo "ERROR: missing required arguments." >&2
    usage >&2
    exit 1
fi

if [[ ! -s "$fasta" ]]; then
    echo "ERROR: FASTA '$fasta' not found or empty." >&2
    exit 1
fi

if [[ ! -s "$gff" ]]; then
    echo "ERROR: GFF '$gff' not found or empty." >&2
    exit 1
fi

mkdir -p "$(dirname "$output")"

tmp_dir="$(mktemp -d)"
trap 'rm -rf "$tmp_dir"' EXIT

bed_file="$tmp_dir/all_blocks.bed"

awk -F'\t' '
    $0 !~ /^#/ {
        block_id = ""
        n = split($9, attrs, ";")
        for (i = 1; i <= n; i++) {
            if (attrs[i] ~ /^ID=/) {
                sub(/^ID=/, "", attrs[i])
                block_id = attrs[i]
            }
        }
        if (block_id != "") {
            start0 = $4 - 1
            end1 = $5
            name = block_id "__" $1 ":" $4 "-" $5
            print $1 "\t" start0 "\t" end1 "\t" name
        }
    }
' "$gff" > "$bed_file"

if [[ ! -s "$bed_file" ]]; then
    echo "ERROR: no intervals were extracted from '$gff'." >&2
    exit 1
fi

bedtools getfasta \
    -fi "$fasta" \
    -bed "$bed_file" \
    -nameOnly \
    -fo "$output"

if [[ ! -s "$output" ]]; then
    echo "ERROR: output FASTA '$output' is empty." >&2
    exit 1
fi
