#!/bin/bash
set -euo pipefail

# Run the preprocessing pipeline for all selected members

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

for member in $(seq 2 5); do
    echo "=============================="
    echo "Processing member ${member}"
    echo "=============================="

    #bash 01_compute_yearly_means.sh "${member}"
    bash 02_compute_ohc_full_depth.sh "${member}"
    bash 03_compute_ohc_layers.sh "${member}"
done

echo "All members processed."
