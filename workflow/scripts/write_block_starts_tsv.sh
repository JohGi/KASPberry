#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  write_block_starts_tsv.sh --input <filtered_blocks.gff> --output <block_starts.tsv>
EOF
}

parse_args() {
  input=""
  output=""

  while [[ $# -gt 0 ]]; do
    case "$1" in
      --input)
        input="${2:-}"
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
        echo "Unknown argument: $1" >&2
        usage >&2
        exit 1
        ;;
    esac
  done

  if [[ -z "$input" || -z "$output" ]]; then
    echo "Both --input and --output are required." >&2
    usage >&2
    exit 1
  fi
}

write_block_starts_tsv() {
  local input_path="$1"
  local output_path="$2"

  mkdir -p "$(dirname "$output_path")"

  {
    printf "block_id\tsample\tblock_start_1based\n"
    awk -F'\t' '!/^#/ && NF {
      split($9, a, "=")
      print a[2], $1, $4
    }' OFS='\t' "$input_path" \
      | sort -t $'\t' -k1,1 -k2,2
  } > "$output_path"
}

main() {
  parse_args "$@"
  write_block_starts_tsv "$input" "$output"
}

main "$@"
