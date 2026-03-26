#!/bin/bash

# ------------------------------------------------------------------
# Assemble historical ensemble OHC with a 3000years spin-up  (members r31–r40)
#
# This script constructs ensemble ocean heat content (OHC) datasets
# from additional historical members (r31–r40) of the CNRM-CM6-1 model.
#
# time segments are first merged into continuous yearly time series (1850–2014),
# then stacked across members to create a four-dimensional dataset:
# (x, y, time, ensemble).
#
# Outputs are produced for:
#   - full-depth OHC
#   - 0–300 m
#   - 300–2000 m
#   - below 2000 m
#
# Final files are written to:
#   /cnrm/ioga/Users/seguy/hist_ens/hist_3000
# ------------------------------------------------------------------


set -euo pipefail
shopt -s nullglob

OUT_DIR="/cnrm/ioga/Users/seguy/hist_ens/hist_3000"
TMP_OUT_DIR="${OUT_DIR}/tmp_members_r31_40"

mkdir -p "${OUT_DIR}"
mkdir -p "${TMP_OUT_DIR}"

PERIODS=(
  "185001-189912"
  "190001-194912"
  "195001-199912"
  "200001-201412"
)

LAYERS=(
  "tot"
  "0-300"
  "300-2000"
  "2000-bottom"
)

echo "======================================"
echo "Building ensemble dataset r31–r40"
echo "Output directory: ${OUT_DIR}"
echo "======================================"

for LAYER in "${LAYERS[@]}"; do

    echo ""
    echo "Processing layer: ${LAYER}"

    case "${LAYER}" in
      tot)
        PREFIX_IN="ohc_2D_yearly"
        OUTFILE="${OUT_DIR}/ohc_2D_hist_3000.nc"
        ;;
      0-300)
        PREFIX_IN="ohc_2D_0-300_yearly"
        OUTFILE="${OUT_DIR}/ohc_2D_hist_3000_0-300.nc"
        ;;
      300-2000)
        PREFIX_IN="ohc_2D_300-2000_yearly"
        OUTFILE="${OUT_DIR}/ohc_2D_hist_3000_300-2000.nc"
        ;;
      2000-bottom)
        PREFIX_IN="ohc_2D_2000-bottom_yearly"
        OUTFILE="${OUT_DIR}/ohc_2D_hist_3000_2000-bottom.nc"
        ;;
    esac

    member_files=()

    for MEMBER in $(seq 31 40); do

        IN_DIR="/cnrm/ioga/Users/seguy/hist_ens/hist_${MEMBER}/tmp"
        MEMBER_OUT="${TMP_OUT_DIR}/${PREFIX_IN}_r${MEMBER}.nc"

        files_to_merge=()

        for PERIOD in "${PERIODS[@]}"; do
            f="${IN_DIR}/${PREFIX_IN}_${PERIOD}_r${MEMBER}.nc"

            if [ ! -f "$f" ]; then
                echo "Missing file: $f"
                exit 1
            fi

            files_to_merge+=("$f")
        done

        echo "Merging member r${MEMBER}"

        cdo -O -f nc4 mergetime \
            "${files_to_merge[@]}" \
            "${MEMBER_OUT}"

        member_files+=("${MEMBER_OUT}")
    done

    echo "Stacking ensemble dimension"

    ncecat -O -4 --mrd -u ensemble \
        "${member_files[@]}" \
        "${OUTFILE}"

    echo "Created:"
    echo "${OUTFILE}"

done

echo ""
echo "======================================"
echo "Done building hist_3000 ensemble files ✅"
echo "======================================"
