#!/usr/bin/env bash
# Wrapper for scripts/build_html.py: pages/full.html + isolated pages/ from content/
#
# Usage:
#   bash scripts/build_site.sh           # writes pages/full.html and pages/*.html
#   bash scripts/build_site.sh --clean   # remove pages/ then rebuild
#
# Requirements:
#   pip install pymdown-extensions pyyaml markdown

set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
cd "$ROOT"

PYTHON="${PYTHON:-python3}"

run_build() {
  echo "--- Building pages/full.html and isolated lesson pages ---"
  "$PYTHON" "$HERE/build_html.py" --root "$ROOT"
}

case "${1:-}" in
  --clean)
    echo "Cleaning pages/..."
    rm -rf pages
    run_build
    ;;
  '')
    run_build
    ;;
  *)
    echo "Usage: bash scripts/build_site.sh [--clean]"
    exit 2
    ;;
esac
