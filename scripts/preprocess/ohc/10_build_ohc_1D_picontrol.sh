#!/bin/bash

set -euo pipefail
shopt -s nullglob

# ------------------------------------------------------------------
# Build global 1D piControl OHC time series from annual 2D OHC files
#
# This script:
#   1. multiplies annual 2D OHC fields (J m^-2) by areacello (m^2),
#   2. sums over x,y to obtain a global 1D OHC time series (J),
#   3. concatenates all annual files along time.
#
# Only the full-column OHC is processed (no depth-layer files).
# ------------------------------------------------------------------

# Resolve repository root
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(dirname "$(dirname "$SCRIPT_DIR")")"

# Metadata
AREACELLO_FILE="${REPO_ROOT}/../data/metadata/areacello_Ofx_CNRM-CM6-1_piControl_r1i1p1f2_gn.nc"

# Input/output directories
IN_DIR="/cnrm/ioga/Users/seguy/pic_ens"
OUT_DIR="/cnrm/ioga/Users/seguyr/pic_ens/pic"
TMP_DIR="${OUT_DIR}/tmp_ohc1D"

mkdir -p "${OUT_DIR}" "${TMP_DIR}"

OUTFILE="${OUT_DIR}/ohc_1D_picontrol.nc"

if [ ! -f "${AREACELLO_FILE}" ]; then
    echo "Missing metadata file: ${AREACELLO_FILE}"
    exit 1
fi

echo "======================================"
echo "Building global 1D piControl OHC"
echo "Input directory : ${IN_DIR}"
echo "Area file       : ${AREACELLO_FILE}"
echo "Output file     : ${OUTFILE}"
echo "======================================"

annual_files=()

for file in "${IN_DIR}"/ohc_2D_picontrol_yearly_*.nc; do
    [ -e "$file" ] || continue

    period=$(basename "$file" | grep -oE '[0-9]{6}-[0-9]{6}')
    tmp_out="${TMP_DIR}/ohc_1D_picontrol_yearly_${period}.nc"

    echo "Processing ${period}"

    # OHC_2D is in J/m²
    # areacello is in m²
    # result after fldsum is in J
    cdo -O -L fldsum -mul "$file" "${AREACELLO_FILE}" "${tmp_out}"

    annual_files+=("${tmp_out}")
done

if [ ${#annual_files[@]} -eq 0 ]; then
    echo "No annual OHC 2D files found in ${IN_DIR}"
    exit 1
fi

echo "Concatenating ${#annual_files[@]} annual files along time"
cdo -O -f nc4 mergetime "${annual_files[@]}" "${OUTFILE}"

echo "======================================"
echo "Done."
echo "Created: ${OUTFILE}"
echo "======================================"
