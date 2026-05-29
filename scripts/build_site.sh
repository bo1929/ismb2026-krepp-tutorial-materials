#!/usr/bin/env bash
# Wrapper for scripts/build_html.py: single-page tutorial.html from content/
#
# Usage:
#   bash scripts/build_site.sh           # writes tutorial.html
#   bash scripts/build_site.sh --clean   # remove tutorial.html then rebuild
#
# Requirements:
#   pip install pymdown-extensions pyyaml markdown

set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
cd "$ROOT"

PYTHON="${PYTHON:-python3}"

run_build() {
  echo "--- Building tutorial.html ---"
  "$PYTHON" "$HERE/build_html.py" --root "$ROOT"
}

case "${1:-}" in
  --clean)
    echo "Cleaning tutorial.html..."
    rm -f tutorial.html
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
