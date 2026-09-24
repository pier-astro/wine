#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

if [[ -z "${HEADAS:-}" ]] || ! command -v initpackage >/dev/null || ! command -v hmake >/dev/null; then
    echo "Initialize HEASoft/XSPEC before building (HEADAS, initpackage, hmake)." >&2
    exit 1
fi

./clean.sh
initpackage wine lmodel.dat .
hmake
echo "Built $PWD/libwine.dylib (or the platform's shared-library suffix)."
