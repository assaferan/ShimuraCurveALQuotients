#!/usr/bin/env bash
# Download the saved weight-3/2 bases (data/dual/README.md) into $DUALW_CACHE (default ~/dualw) and check each
# against data/dual/bases.md5.  Usage: tools/fetch_dual_bases.sh [basis_M_B.txt ...]   (default: all of them)
set -euo pipefail
ZENODO_URL="${ZENODO_URL:-https://zenodo.org/records/RECORD_ID/files}"   # set once the deposit exists
DEST="${DUALW_CACHE:-$HOME/dualw}"
LIST="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/data/dual/bases.md5"
mkdir -p "$DEST"
md5of() { if command -v md5sum >/dev/null; then md5sum "$1" | cut -c1-32; else md5 -q "$1"; fi; }
want=("$@")
grep -v '^#' "$LIST" | while read -r sum file size; do
    [ -z "$file" ] && continue
    if [ ${#want[@]} -gt 0 ] && [[ ! " ${want[*]} " =~ " $file " ]]; then continue; fi
    if [ -f "$DEST/$file" ] && [ "$(md5of "$DEST/$file")" = "$sum" ]; then echo "have  $file"; continue; fi
    echo "fetch $file ($size bytes)"
    curl -fL --retry 3 -o "$DEST/$file.part" "$ZENODO_URL/$file?download=1"
    got="$(md5of "$DEST/$file.part")"
    if [ "$got" != "$sum" ]; then echo "checksum mismatch for $file: $got, expected $sum" >&2; rm -f "$DEST/$file.part"; exit 1; fi
    mv "$DEST/$file.part" "$DEST/$file"
done
