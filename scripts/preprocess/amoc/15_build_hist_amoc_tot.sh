#!/bin/bash


set -euo pipefail
shopt -s nullglob

# ------------------------------------------------------------------
# Assemble historical AMOC members r2-r30 into one ensemble file
# ------------------------------------------------------------------

IN_DIR="/cnrm/ioga/Users/seguy/hist_ens/amoc_members"
OUT_DIR="/cnrm/ioga/Users/seguy/hist_ens/hist_tot"

mkdir -p "${OUT_DIR}"

OUTFILE="${OUT_DIR}/amoc_hist_tot.nc"

member_files=()

for MEMBER in $(seq 2 5); do
    f="${IN_DIR}/amoc_hist_r${MEMBER}.nc"

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
