#!/bin/bash
#
set -euo pipefail
shopt -s nullglob

# ------------------------------------------------------------------
# Assemble piControl AMOC pseudo-members r1-r29 into one ensemble file
# ------------------------------------------------------------------

IN_DIR="/cnrm/ioga/Users/seguyr/pic_ens/pic_amoc"
OUT_DIR="/cnrm/ioga/Users/seguyr/pic_ens/pic_tot"

mkdir -p "${OUT_DIR}"

OUTFILE="${OUT_DIR}/amoc_pic_tot.nc"

member_files=()

for MEMBER in $(seq 1 5); do
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
