#!/usr/bin/env bash
set -euo pipefail

# Filter SibeliaZ GFF blocks using a two-pass awk strategy.
# Keep blocks where:
#   - every feature length is greater than or equal to min_len
#   - number of distinct sequence names is at least nb_samples
#   - no sequence name appears more than once in the same block
#   - all features in the block have the same strand
#
# Output: filtered GFF written to the output file.

usage() {
  cat <<'EOF'
Usage:
  filter_sibeliaz_blocks.sh --gff blocks.gff --nb_samples N [--min_len 500] --output filtered.gff

Arguments:
  --gff         Input GFF file
  --nb_samples Required minimum number of distinct sequence names per block
  --min_len    Minimum feature length required for every block occurrence, inclusive (default: 500)
  --output     Output filtered GFF file
EOF
}

parse_args() {
  gff=""
  nb_samples=""
  min_len=500
  output=""

  while [[ $# -gt 0 ]]; do
    case "$1" in
      --gff)
        gff="${2:-}"
        shift 2
        ;;
      --nb_samples)
        nb_samples="${2:-}"
        shift 2
        ;;
      --min_len)
        min_len="${2:-}"
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
        echo "ERROR: Unknown argument: $1" >&2
        usage >&2
        exit 2
        ;;
    esac
  done

  if [[ -z "$gff" ]]; then
    echo "ERROR: --gff is required." >&2
    exit 2
  fi

  if [[ -z "$nb_samples" ]]; then
    echo "ERROR: --nb_samples is required." >&2
    exit 2
  fi

  if [[ -z "$output" ]]; then
    echo "ERROR: --output is required." >&2
    exit 2
  fi

  if [[ ! -f "$gff" ]]; then
    echo "ERROR: GFF not found: $gff" >&2
    exit 2
  fi

  if ! [[ "$nb_samples" =~ ^[0-9]+$ ]]; then
    echo "ERROR: --nb_samples must be a non-negative integer." >&2
    exit 2
  fi

  if ! [[ "$min_len" =~ ^[0-9]+$ ]]; then
    echo "ERROR: --min_len must be a non-negative integer." >&2
    exit 2
  fi
}

run_filter() {
  local gff="$1"
  local min_len="$2"
  local nb_samples="$3"
  local output="$4"

  mkdir -p "$(dirname "$output")"

  awk -F'\t' -v min_len="$min_len" -v nb_samples="$nb_samples" '
    FNR == NR {
      if (!match($9, /(^|;)ID=([^;]+)/, m)) {
        next
      }

      id = m[2]
      sequence_name = $1
      length_bp = $5 - $4 + 1

      if (!(id in min_feature_len) || length_bp < min_feature_len[id]) {
        min_feature_len[id] = length_bp
      }

      key = id SUBSEP sequence_name
      feature_count[key]++

      if (feature_count[key] == 1) {
        n_sequences[id]++
      }

      if (feature_count[key] > 1) {
        has_paralog[id] = 1
      }

      strand = $7
      if (!(id in first_strand)) {
        first_strand[id] = strand
      } else if (first_strand[id] != strand) {
        mixed_strand[id] = 1
      }

      next
    }

    {
      if (!match($9, /(^|;)ID=([^;]+)/, m)) {
        next
      }

      id = m[2]

      if (min_feature_len[id] < min_len) {
        next
      }

      if (n_sequences[id] < nb_samples) {
        next
      }

      if (has_paralog[id]) {
        next
      }

      if (mixed_strand[id]) {
        next
      }

      print
    }
  ' "$gff" "$gff" > "$output"
}

main() {
  parse_args "$@"
  run_filter "$gff" "$min_len" "$nb_samples" "$output"
}

main "$@"