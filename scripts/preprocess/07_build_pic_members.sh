#!/bin/bash


set -euo pipefail
shopt -s nullglob

# ------------------------------------------------------------------
# Build 39 piControl pseudo-members from branching years
#
# This script extracts 39 piControl segments of length 165 years from
# yearly OHC diagnostics previously computed from the CNRM-CM6-1
# piControl simulation.
#
# The start years are read from a NetCDF file containing a variable
# named "year". For each start year, a 165-year segment is extracted,
# concatenated in time, and reindexed onto the historical time axis
# 1850-2014.
#
# Outputs are produced for:
#   - 3D OHC
#   - full-column 2D OHC
#   - 0-300 m OHC
#   - 300-2000 m OHC
#   - 2000-bottom OHC
# ------------------------------------------------------------------


SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(dirname "$(dirname "$SCRIPT_DIR")")"


BRANCH_FILE="${REPO_ROOT}/data/metadata/hist_start_years.nc"


PIC_INPUT_DIR="/cnrm/ioga/Users/seguy/pic_ens"
OUT_DIR="/cnrm/ioga/Users/seguy/pic_ens/pic_members"
TMP_DIR="${OUT_DIR}/tmp"

mkdir -p "${OUT_DIR}" "${TMP_DIR}"

# Read branching years from NetCDF variable "year"
mapfile -t START_YEARS < <(
python3 - <<PY
import xarray as xr
ds = xr.open_dataset("${BRANCH_FILE}")
years = ds["year"].values.tolist()
for y in years:
    print(int(y))
PY
)

N_MEMBERS="${#START_YEARS[@]}"

echo "======================================"
echo "Building piControl pseudo-members"
echo "Branching-year file : ${BRANCH_FILE}"
echo "Input directory     : ${PIC_INPUT_DIR}"
echo "Output directory    : ${OUT_DIR}"
echo "Number of members   : ${N_MEMBERS}"
echo "======================================"

# Input/output prefixes
IN_PREFIXES=(
  "ohc_3D_picontrol_yearly"
  "ohc_2D_picontrol_yearly"
  "ohc_2D_0-300_picontrol_yearly"
  "ohc_2D_300-2000_picontrol_yearly"
  "ohc_2D_2000-bottom_picontrol_yearly"
)

OUT_PREFIXES=(
  "ohc_3D_pic_yearly"
  "ohc_2D_pic_yearly"
  "ohc_2D_0-300_pic_yearly"
  "ohc_2D_300-2000_pic_yearly"
  "ohc_2D_2000-bottom_pic_yearly"
)

# Loop over pseudo-members
for idx in "${!START_YEARS[@]}"; do
    MEMBER_ID=$((idx + 1))
    START_YEAR="${START_YEARS[$idx]}"
    END_YEAR=$((START_YEAR + 164))

    echo "--------------------------------------"
    echo "Building pic_r${MEMBER_ID}"
    echo "Source years : ${START_YEAR}-${END_YEAR}"
    echo "Target years : 1850-2014"
    echo "--------------------------------------"

    for p in "${!IN_PREFIXES[@]}"; do
        IN_PREFIX="${IN_PREFIXES[$p]}"
        OUT_PREFIX="${OUT_PREFIXES[$p]}"

        files_to_merge=()

        for year in $(seq "${START_YEAR}" "${END_YEAR}"); do
            period=$(printf "%04d01-%04d12" "${year}" "${year}")
            f="${PIC_INPUT_DIR}/${IN_PREFIX}_${period}.nc"

            if [ ! -f "$f" ]; then
                echo "Missing file for pic_r${MEMBER_ID}: $f"
                exit 1
            fi

            files_to_merge+=("$f")
        done

        merged_file="${TMP_DIR}/${OUT_PREFIX}_r${MEMBER_ID}_raw.nc"
        outfile="${OUT_DIR}/${OUT_PREFIX}_r${MEMBER_ID}.nc"

        echo "Merging ${OUT_PREFIX} for pic_r${MEMBER_ID}"

        cdo -O -f nc4 mergetime "${files_to_merge[@]}" "${merged_file}"

        # Reindex time axis to 1850-2014 with 1-year increment
        cdo -O -f nc4 settaxis,1850-01-01,00:00:00,1year "${merged_file}" "${outfile}"

        rm -f "${merged_file}"
    done
done

echo "======================================"
echo "All piControl pseudo-members created"
echo "Outputs written to: ${OUT_DIR}"
echo "======================================"
