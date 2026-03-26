#!/bin/bash

set -euo pipefail
shopt -s nullglob

# ------------------------------------------------------------------
# Assemble piControl pseudo-members r2-r30 into ensemble datasets
#
# This script stacks piControl pseudo-members (r2 to r30) along a new
# ensemble dimension to produce four-dimensional NetCDF datasets:
# (x, y, time, ensemble).
#
# Outputs are generated for:
#   - full-depth OHC
#   - 0-300 m
#   - 300-2000 m
#   - below 2000 m
#
# Final files are written to:
#   /cnrm/ioga/Users/seguy/CMIP6/PiControl/pic_tot
# ------------------------------------------------------------------

IN_DIR="/cnrm/ioga/Users/seguy/pic_ens/pic_members"
OUT_DIR="/cnrm/ioga/Users/seguy/pic_ens/pic_tot"

mkdir -p "${OUT_DIR}"

LAYERS=(
  "tot"
  "0-300"
  "300-2000"
  "2000-bottom"
)

echo "======================================"
echo "Building piControl ensemble r2-r30"
echo "Input directory : ${IN_DIR}"
echo "Output directory: ${OUT_DIR}"
echo "======================================"

for LAYER in "${LAYERS[@]}"; do

    case "${LAYER}" in
      tot)
        PREFIX_IN="ohc_2D_pic_yearly"
        OUTFILE="${OUT_DIR}/ohc_2D_pic_tot.nc"
        ;;
      0-300)
        PREFIX_IN="ohc_2D_0-300_pic_yearly"
        OUTFILE="${OUT_DIR}/ohc_2D_pic_tot_0-300.nc"
        ;;
      300-2000)
        PREFIX_IN="ohc_2D_300-2000_pic_yearly"
        OUTFILE="${OUT_DIR}/ohc_2D_pic_tot_300-2000.nc"
        ;;
      2000-bottom)
        PREFIX_IN="ohc_2D_2000-bottom_pic_yearly"
        OUTFILE="${OUT_DIR}/ohc_2D_pic_tot_2000-bottom.nc"
        ;;
    esac

    echo "--------------------------------------"
    echo "Processing layer: ${LAYER}"
    echo "--------------------------------------"

    member_files=()

    for MEMBER in $(seq 2 30); do
        f="${IN_DIR}/${PREFIX_IN}_r${MEMBER}.nc"

        if [ ! -f "$f" ]; then
            echo "Missing file: $f"
            exit 1
        fi

        member_files+=("$f")
    done

    ncecat -O -4 --mrd -u ensemble \
        "${member_files[@]}" \
        "${OUTFILE}"

    echo "Created: ${OUTFILE}"
done

echo "======================================"
echo "Done building pic_tot files"
echo "======================================"
