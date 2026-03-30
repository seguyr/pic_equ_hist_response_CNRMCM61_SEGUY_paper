#!/bin/bash
set -euo pipefail
shopt -s nullglob

# ------------------------------------------------------------------
# Compute yearly ocean heat content from CNRM-CM6-1 piControl output
#
# This script computes yearly ocean heat content (OHC) diagnostics
# from CNRM-CM6-1 piControl thetao and thkcello outputs.
#
# Because thetao and thkcello chunking is not temporally consistent,
# the workflow:
#   1. computes yearly means,
#   2. splits them into one file per year,
#   3. matches thetao and thkcello year by year,
#   4. computes 3D and vertically integrated OHC.
#
# Output files preserve the original simulation year, e.g.:
#   185001-185012, 185101-185112, ...
# ------------------------------------------------------------------

THETAO_SRC="/cnrm/cmip6/CMIP/CNRM-CERFACS/CNRM-CM6-1/piControl/r1i1p1f2/Omon/thetao/gn/latest"
THKCELLO_SRC="/cnrm/cmip6/CMIP/CNRM-CERFACS/CNRM-CM6-1/piControl/r1i1p1f2/Omon/thkcello/gn/latest"

OUT_DIR="/cnrm/ioga/Users/seguy/pic_ens"
YEARLY_DIR="${OUT_DIR}/yearly"
YEAR_SPLIT_DIR="${OUT_DIR}/year_split"
TMP_DIR="${OUT_DIR}/tmp"

THETAO_YEAR_DIR="${YEAR_SPLIT_DIR}/thetao"
THKCELLO_YEAR_DIR="${YEAR_SPLIT_DIR}/thkcello"

mkdir -p "${OUT_DIR}" "${YEARLY_DIR}" "${YEAR_SPLIT_DIR}" "${TMP_DIR}"
mkdir -p "${THETAO_YEAR_DIR}" "${THKCELLO_YEAR_DIR}"

# Physical constants
CP=3992
RHO=1025
CPRHO=$(awk "BEGIN {print ${CP}*${RHO}}")

echo "======================================"
echo "Computing piControl OHC"
echo "thetao source   : ${THETAO_SRC}"
echo "thkcello source : ${THKCELLO_SRC}"
echo "output dir      : ${OUT_DIR}"
echo "Cp*rho          : ${CPRHO}"
echo "======================================"

# --------------------------------------------------
# Step 1: compute yearly means for thetao
# --------------------------------------------------
cd "${THETAO_SRC}" || exit 1

for filet in thetao_*.nc; do
    [ -e "$filet" ] || continue

    clean_file="${filet//@/}"
    yearly_file="${YEARLY_DIR}/${clean_file/Omon/yearly}"

    if [ ! -f "$yearly_file" ]; then
        echo "Yearly mean thetao: $filet"
        cdo -L yearmean "$filet" "$yearly_file"
    else
        echo "Already exists, skipping: $yearly_file"
    fi
done

# --------------------------------------------------
# Step 2: compute yearly means for thkcello
# --------------------------------------------------
cd "${THKCELLO_SRC}" || exit 1

for filet in thkcello_*.nc; do
    [ -e "$filet" ] || continue

    clean_file="${filet//@/}"
    yearly_file="${YEARLY_DIR}/${clean_file/Omon/yearly}"

    if [ ! -f "$yearly_file" ]; then
        echo "Yearly mean thkcello: $filet"
        cdo -L yearmean "$filet" "$yearly_file"
    else
        echo "Already exists, skipping: $yearly_file"
    fi
done

# --------------------------------------------------
# Step 3: split yearly files into one file per year
# --------------------------------------------------
cd "${YEARLY_DIR}" || exit 1

for file in thetao_yearly_*.nc; do
    [ -e "$file" ] || continue
    echo "Splitting thetao yearly file: $file"
    cdo splityear "$file" "${THETAO_YEAR_DIR}/thetao_year_"
done

for file in thkcello_yearly_*.nc; do
    [ -e "$file" ] || continue
    echo "Splitting thkcello yearly file: $file"
    cdo splityear "$file" "${THKCELLO_YEAR_DIR}/thkcello_year_"
done

# --------------------------------------------------
# Step 4: rename split files using the real simulation year
# --------------------------------------------------
for f in "${THETAO_YEAR_DIR}"/thetao_year_*.nc; do
    [ -e "$f" ] || continue
    y=$(cdo -s showyear "$f" | awk '{print $1}')
    target="${THETAO_YEAR_DIR}/thetao_year_${y}.nc"
    if [ "$f" != "$target" ]; then
        mv -f "$f" "$target"
    fi
done

for f in "${THKCELLO_YEAR_DIR}"/thkcello_year_*.nc; do
    [ -e "$f" ] || continue
    y=$(cdo -s showyear "$f" | awk '{print $1}')
    target="${THKCELLO_YEAR_DIR}/thkcello_year_${y}.nc"
    if [ "$f" != "$target" ]; then
        mv -f "$f" "$target"
    fi
done

# --------------------------------------------------
# Step 5: compute OHC year by year
# --------------------------------------------------
cd "${THETAO_YEAR_DIR}" || exit 1

for thetao_file in thetao_year_*.nc; do
    [ -e "$thetao_file" ] || continue

    year=$(cdo -s showyear "$thetao_file" | awk '{print $1}')
    thkcello_file="${THKCELLO_YEAR_DIR}/thkcello_year_${year}.nc"

    if [ ! -f "$thkcello_file" ]; then
        echo "Skipping year ${year}: missing thkcello file"
        continue
    fi

    period="${year}01-${year}12"
    thetaoK_file="${TMP_DIR}/thetaoK_year_${year}.nc"

    ohc3D_file="${OUT_DIR}/ohc_3D_picontrol_yearly_${period}.nc"
    ohc2D_file="${OUT_DIR}/ohc_2D_picontrol_yearly_${period}.nc"
    ohc2D_0300_file="${OUT_DIR}/ohc_2D_0-300_picontrol_yearly_${period}.nc"
    ohc2D_3002000_file="${OUT_DIR}/ohc_2D_300-2000_picontrol_yearly_${period}.nc"
    ohc2D_2000bottom_file="${OUT_DIR}/ohc_2D_2000-bottom_picontrol_yearly_${period}.nc"

    echo "--------------------------------------"
    echo "Processing year   : ${year}"
    echo "thetao file       : ${thetao_file}"
    echo "thkcello file     : ${thkcello_file}"
    echo "output period     : ${period}"
    echo "--------------------------------------"

    # Convert temperature from °C to K
    cdo addc,273.15 "${thetao_file}" "${thetaoK_file}"

    # 3D OHC in J/m²
    cdo mulc,"$CPRHO" -mul "${thetaoK_file}" "${thkcello_file}" "${ohc3D_file}"

    # Full-column 2D OHC in J/m²
    cdo vertsum "${ohc3D_file}" "${ohc2D_file}"

    # Layer-integrated 2D OHC in J/m²
    cdo vertsum -sellevidx,1/34  "${ohc3D_file}" "${ohc2D_0300_file}"
    cdo vertsum -sellevidx,34/54 "${ohc3D_file}" "${ohc2D_3002000_file}"
    cdo vertsum -sellevidx,55/75 "${ohc3D_file}" "${ohc2D_2000bottom_file}"

    rm -f "${thetaoK_file}"
done

echo "======================================"
echo "piControl OHC computation completed"
echo "Outputs written to: ${OUT_DIR}"
echo "======================================"
