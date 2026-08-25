#!/usr/bin/env bash
set -euo pipefail

MFEPRIMER_VERSION="4.5.1"
MFEPRIMER_URL="https://github.com/quwubin/MFEprimer-3.0/releases/download/v4.5.1/mfeprimer-4.5.1-linux-amd64.gz"
MFEPRIMER_SHA256="56abb0789497a6273e0b7d226671d5f1d959d3bceaba5be53f0d2b2edc13c839"

archive="$(mktemp)"
binary_tmp="$(mktemp "${CONDA_PREFIX}/bin/.mfeprimer.XXXXXX")"

cleanup() {
    rm -f "$archive" "$binary_tmp"
}
trap cleanup EXIT

curl --fail --location --retry 3 --output "$archive" "$MFEPRIMER_URL"

if ! printf '%s  %s\n' "$MFEPRIMER_SHA256" "$archive" | sha256sum --check --status; then
    printf 'MFEprimer %s checksum verification failed.\n' "$MFEPRIMER_VERSION" >&2
    exit 1
fi

gzip --decompress --stdout "$archive" > "$binary_tmp"
chmod +x "$binary_tmp"
"$binary_tmp" --help >/dev/null
mv -f "$binary_tmp" "$CONDA_PREFIX/bin/mfeprimer"
