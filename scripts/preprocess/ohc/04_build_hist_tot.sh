#!/bin/bash
# ------------------------------------------------------------------
# Assemble historical ensemble OHC with less than 1000 years of spin-up
# with all historical CMIP6 members (members r2–r30)
# r1 thetao variable was not complete that is why we begin at number 2
#
# This script concatenates vertically integrated ocean heat content
# (OHC) fields from ensemble members r2–r30 of the CNRM-CM6-1
# historical simulations.
#
# For each member, time segments (1850–1899, 1900–1949,
# 1950–1999, 2000–2014) are merged to produce a continuous
# yearly time series (1850–2014; 165 years).
#
# The resulting member time series are then stacked along a new
# ensemble dimension using NCO to produce a four-dimensional dataset:
# (x, y, time, ensemble).
#
# Outputs are generated for:
#   - full-depth OHC
#   - 0–300 m
#   - 300–2000 m
#   - below 2000 m
#
# Final files are written to:
#   /cnrm/ioga/Users/seguy/hist_ens/hist_tot
# ------------------------------------------------------------------
#
#
#
set -euo pipefail
shopt -s nullglob

OUT_DIR="/cnrm/ioga/Users/seguy/hist_ens/hist_tot"
TMP_OUT_DIR="${OUT_DIR}/tmp_members"

mkdir -p "${OUT_DIR}"
mkdir -p "${TMP_OUT_DIR}"

# Time chunks
PERIODS=(
  "185001-189912"
  "190001-194912"
  "195001-199912"
  "200001-201412"
)

# Layers to process
LAYERS=(
  "tot"
  "0-300"
  "300-2000"
  "2000-bottom"
)

echo "======================================"
echo "Building ensemble historical dataset"
echo "Members: r2 to r30"
echo "Output dir: ${OUT_DIR}"
echo "======================================"

for LAYER in "${LAYERS[@]}"; do

    echo ""
    echo "Processing layer: ${LAYER}"
    echo "--------------------------------------"

    case "${LAYER}" in
      tot)
        PREFIX_IN="ohc_2D_yearly"
        OUTFILE="${OUT_DIR}/ohc_2D_hist_tot.nc"
        ;;
      0-300)
        PREFIX_IN="ohc_2D_0-300_yearly"
        OUTFILE="${OUT_DIR}/ohc_2D_hist_tot_0-300.nc"
        ;;
      300-2000)
        PREFIX_IN="ohc_2D_300-2000_yearly"
        OUTFILE="${OUT_DIR}/ohc_2D_hist_tot_300-2000.nc"
        ;;
      2000-bottom)
        PREFIX_IN="ohc_2D_2000-bottom_yearly"
        OUTFILE="${OUT_DIR}/ohc_2D_hist_tot_2000-bottom.nc"
        ;;
    esac

    member_files=()

    # Step 1: merge time for each member
    for MEMBER in $(seq 2 5); do

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


    # Step 2: merge ensemble dimension
    echo "Stacking ensemble dimension"

    ncecat -O -4 --mrd -u ensemble \
        "${member_files[@]}" \
        "${OUTFILE}"

    echo "Created:"
    echo "${OUTFILE}"

done

echo ""
echo "======================================"
echo "All layers processed successfully ✅"
echo "======================================"
