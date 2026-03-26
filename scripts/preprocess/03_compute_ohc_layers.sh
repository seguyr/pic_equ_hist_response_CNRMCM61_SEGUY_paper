#!/bin/bash

# ------------------------------------------------------------------
# Compute layer-integrated ocean heat content
#
# This script computes vertically integrated ocean heat content
# (OHC) within selected depth ranges from previously generated
# 3D OHC fields.
#
# Depth-integrated diagnostics are produced for:
#   - 0–300 m
#   - 300–2000 m
#   - below 2000 m
#
# Output files retain the original temporal segmentation
# (1850–1899, 1900–1949, 1950–1999, 2000–2014) and are saved
# separately for each ensemble member.
# ------------------------------------------------------------------


set -euo pipefail

# Compute layered ocean heat content for one historical member
# Example:
# bash scripts/preprocess/03_compute_ohc_layers.sh 25

MEMBER="${1:?Provide member number (example: 1, 12, 25)}"

TMP_DIR="/cnrm/ioga/Users/seguy/hist_ens/hist_${MEMBER}/tmp"

cd "${TMP_DIR}" || exit 1

echo "======================================"
echo "Computing vertical OHC layers for member r${MEMBER}"
echo "Working directory: ${TMP_DIR}"
echo "======================================"

for file in ohc_3D_yearly_*_r${MEMBER}.nc; do
    [ -e "$file" ] || continue

    period=$(basename "$file" | grep -oE '[0-9]{6}-[0-9]{6}')

    echo "Processing period ${period}, member r${MEMBER}"
    echo "Input file: $file"

    cdo vertsum -sellevidx,1/34 \
        "$file" \
        "ohc_2D_0-300_yearly_${period}_r${MEMBER}.nc"

    cdo vertsum -sellevidx,34/54 \
        "$file" \
        "ohc_2D_300-2000_yearly_${period}_r${MEMBER}.nc"

    cdo vertsum -sellevidx,55/75 \
        "$file" \
        "ohc_2D_2000-bottom_yearly_${period}_r${MEMBER}.nc"
done

echo "Layer-integrated OHC completed for member r${MEMBER}"


