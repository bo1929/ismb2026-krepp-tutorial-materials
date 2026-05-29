#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")/.."

wget --no-check-certificate https://ter-trees.ucsd.edu/data/krepp/index-WoLv1-tiny.tar
mkdir -p data
tar -xf index-WoLv1-tiny.tar -C data/
mv index-WoLv1-tiny.tar data/

mkdir -p results
