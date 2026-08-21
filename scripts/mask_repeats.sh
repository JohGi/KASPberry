#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Mask repeats in a block FASTA using RepeatMasker.

Usage:
  mask_repeats.sh \
    --fasta input.fasta \
    [--te-lib repeats.fasta] \
    --outdir masked_dir \
    --threads 3 \
    --output output.masked
EOF
}

fasta=""
te_lib=""
outdir=""
threads=""
output=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --fasta)
            fasta="${2:-}"
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
        --output)
            output="${2:-}"
            shift 2
            ;;
        *)
            echo "Error: unknown argument '$1'." >&2
            usage >&2
            exit 1
            ;;
    esac
done

if [[ -z "$fasta" || -z "$outdir" || -z "$threads" || -z "$output" ]]; then
    echo "Error: missing required arguments." >&2
    usage >&2
    exit 1
fi

if [[ ! -s "$fasta" ]]; then
    echo "Error: FASTA '$fasta' not found or empty." >&2
    exit 1
fi

if [[ -n "$te_lib" && ! -s "$te_lib" ]]; then
    echo "Error: TE library '$te_lib' not found or empty." >&2
    exit 1
fi

mkdir -p "$outdir"
mkdir -p "$(dirname "$output")"

block_name="$(basename "$fasta")"
block_tmp_dir="$(mktemp -d "${outdir}/repeatmasker_${block_name}.XXXXXX")"
trap 'rm -rf "$block_tmp_dir"' EXIT

masked_candidate="$block_tmp_dir/${block_name}.masked"

repeatmasker_args=(
    -pa "$threads"
    -no_is
    -dir "$block_tmp_dir"
)

if [[ -n "$te_lib" ]]; then
    repeatmasker_args+=( -lib "$te_lib" )
else
    # RepeatMasker >= 4.2 requires either a custom library or a configured
    # FamDB even when -noint is used. KASPberry therefore supplies a small
    # technical placeholder library solely to satisfy this requirement.
    # -noint ensures that only simple repeats and low-complexity regions
    # are actually masked.
    script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    placeholder_lib="${script_dir}/resources/repeatmasker_placeholder.fa"

    if [[ ! -s "$placeholder_lib" ]]; then
        echo "Error: RepeatMasker placeholder library '$placeholder_lib' not found or empty." >&2
        exit 1
    fi

    repeatmasker_args+=(
        -lib "$placeholder_lib"
        -noint
    )
fi

RepeatMasker "${repeatmasker_args[@]}" "$fasta"

if [[ -f "$masked_candidate" ]]; then
    cp "$masked_candidate" "$output"
else
    cp "$fasta" "$output"
fi
