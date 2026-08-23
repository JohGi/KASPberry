#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  detect_snps.sh \
    --align-dir DIR \
    --output FILE \
    --min-flank INT

Description:
  Prepare temporary alignment files for seqtui by replacing N/n with '-'
  in sequence lines only, then run seqtui on the temporary files.

Arguments:
  --align-dir   Directory containing input alignment FASTA files
  --output      Output VCF path
  --min-flank   Minimum flank value passed to seqtui
EOF
}

parse_args() {
  align_dir=""
  output=""
  min_flank=""

  while [[ $# -gt 0 ]]; do
    case "$1" in
      --align-dir)
        align_dir="${2:-}"
        shift 2
        ;;
      --output)
        output="${2:-}"
        shift 2
        ;;
      --min-flank)
        min_flank="${2:-}"
        shift 2
        ;;
      --help|-h)
        usage
        exit 0
        ;;
      *)
        echo "Error: unknown argument: $1" >&2
        usage >&2
        exit 1
        ;;
    esac
  done

  if [[ -z "$align_dir" || -z "$output" || -z "$min_flank" ]]; then
    echo "Error: missing required arguments." >&2
    usage >&2
    exit 1
  fi

  if [[ ! -d "$align_dir" ]]; then
    echo "Error: alignment directory not found: $align_dir" >&2
    exit 1
  fi
}

cleanup_tmp_dir() {
  if [[ -n "${tmp_dir:-}" && -d "${tmp_dir:-}" ]]; then
    rm -rf "$tmp_dir"
  fi
}

convert_alignment_for_seqtui() {
  local input_fasta="$1"
  local output_fasta="$2"

  perl -pe 'if (!/^>/) { s/[nN]/-/g }' "$input_fasta" > "$output_fasta"
}

main() {
  parse_args "$@"

  tmp_dir=""
  trap cleanup_tmp_dir EXIT

  tmp_dir="$(mktemp -d)"

  shopt -s nullglob
  local input_files=("$align_dir"/*.aln.fasta)
  shopt -u nullglob

  if [[ ${#input_files[@]} -eq 0 ]]; then
    echo "Error: no alignment files matching *.aln.fasta found in $align_dir" >&2
    exit 1
  fi

  local tmp_files=()
  local input_file=""
  local tmp_file=""

  for input_file in "${input_files[@]}"; do
    tmp_file="$tmp_dir/$(basename "$input_file")"
    convert_alignment_for_seqtui "$input_file" "$tmp_file"
    tmp_files+=("$tmp_file")
  done

  seqtui \
    -v "$min_flank" \
    -d ":" \
    -o "$output" \
    "${tmp_files[@]}"
}

main "$@"
