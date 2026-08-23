#!/usr/bin/env bash
set -euo pipefail

# Rename the header of a single-sequence FASTA and write the result to an output file.

usage() {
  cat <<'EOF'
Usage:
  rename_fasta_header.sh --fasta INPUT.fasta --name NAME --output OUTPUT.fasta

Arguments:
  --fasta    Input FASTA file (must be a single-sequence FASTA)
  --name     New sequence name
  --output   Output FASTA file
EOF
}

parse_args() {
  fasta=""
  name=""
  output=""

  while [[ $# -gt 0 ]]; do
    case "$1" in
      --fasta) fasta="${2:-}"; shift 2 ;;
      --name) name="${2:-}"; shift 2 ;;
      --output) output="${2:-}"; shift 2 ;;
      -h|--help) usage; exit 0 ;;
      *) echo "ERROR: Unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
  done

  if [[ -z "$fasta" || -z "$name" || -z "$output" ]]; then
    echo "ERROR: --fasta, --name and --output are required." >&2
    usage >&2
    exit 1
  fi

  if [[ ! -f "$fasta" ]]; then
    echo "ERROR: FASTA not found: $fasta" >&2
    exit 1
  fi

  if [[ "$name" == '>'* ]]; then
    echo "ERROR: --name must not start with '>'." >&2
    exit 1
  fi

  if [[ "$(grep -c '^>' "$fasta")" -ne 1 ]]; then
    echo "ERROR: FASTA must contain exactly one sequence: $fasta" >&2
    exit 1
  fi
}

rename_fasta_header() {
  local fasta="$1"
  local name="$2"
  local output="$3"

  mkdir -p "$(dirname "$output")"

  awk -v header="$name" '
    BEGIN { done=0 }
    /^>/ && !done { print ">" header; done=1; next }
    { print }
  ' "$fasta" > "$output"
}

main() {
  parse_args "$@"
  rename_fasta_header "$fasta" "$name" "$output"
}

main "$@"
