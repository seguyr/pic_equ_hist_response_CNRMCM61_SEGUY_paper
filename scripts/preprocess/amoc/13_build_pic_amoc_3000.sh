#!/bin/bash

set -euo pipefail
shopt -s nullglob

# ------------------------------------------------------------------
# Assemble piControl AMOC pseudo-members r30-r39 into one ensemble file
# ------------------------------------------------------------------

IN_DIR="/cnrm/ioga/Users/seguy/pic_ens/pic_amoc"
OUT_DIR="/cnrm/ioga/Users/seguy/pic_ens/pic_3000"

mkdir -p "${OUT_DIR}"

OUTFILE="${OUT_DIR}/amoc_pic_3000.nc"

member_files=()

for MEMBER in $(seq 30 39); do
    f="${IN_DIR}/amoc_pic_r${MEMBER}.nc"

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
