#!/bin/bash

# ------------------------------------------------------------------
# Compute vertically integrated ocean heat content (OHC)
#
# This script computes ocean heat content (OHC) from yearly mean
# temperature (thetao) and layer thickness (thkcello) fields for
# each historical ensemble member of the CNRM-CM6-1 model.
#
# Temperature fields are converted from °C to Kelvin and multiplied
# by layer thickness and the volumetric heat capacity of seawater
# (rho * Cp) to obtain 3D OHC fields (J m^-2).
#
# Vertically integrated 2D OHC fields are also produced over the
# full water column for the periods:
#   1850–1899, 1900–1949, 1950–1999, 2000–2014
#
# Outputs are stored per ensemble member and per time chunk.
# ------------------------------------------------------------------


set -euo pipefail

# Compute ocean heat content in J/m² for one historical member.
#
# Outputs:
# - heatc3D_* : 3D fields in J/m²
# - heatc2D_* : 2D vertically integrated fields in J/m² over the full water column
#
# Usage:
#   bash scripts/preprocess/02_compute_ohc_full_depth.sh 1
#   bash scripts/preprocess/02_compute_ohc_full_depth.sh 25
#
#
MEMBER="${1:?Provide member number (example: 1, 12, 25)}"

BASE_DIR="/cnrm/ioga/Users/seguy/hist_ens/hist_${MEMBER}"
TMP_DIR="${BASE_DIR}/tmp"

# Physical constants
CP=3992
RHO=1025
CPRHO=$(awk "BEGIN {print ${CP}*${RHO}}")

mkdir -p "${TMP_DIR}"

cd "${BASE_DIR}" || exit 1

echo "======================================"
echo "Computing OHC for member r${MEMBER}"
echo "Base directory: ${BASE_DIR}"
echo "Temporary directory: ${TMP_DIR}"
echo "Cp*rho = ${CPRHO}"
echo "======================================"

# --------------------------------------------------
# Step 1: merge thetao 25-year files into 50-year files
# to match thkcello temporal chunks
# --------------------------------------------------

echo "Merging thetao 25-year files into 50-year files"

# Expected historical chunks:
# 1850-1874 + 1875-1899 -> 1850-1899
# 1900-1924 + 1925-1949 -> 1900-1949
# 1950-1974 + 1975-1999 -> 1950-1999
# 2000-2014 stays unchanged

if [ -f thetao_yearly_CNRM-CM6-1_historical_r${MEMBER}i1p1f2_gn_185001-187412.nc ] && [ -f thetao_yearly_CNRM-CM6-1_historical_r${MEMBER}i1p1f2_gn_187501-189912.nc ]; then
    cdo mergetime thetao_yearly_CNRM-CM6-1_historical_r${MEMBER}i1p1f2_gn_185001-187412.nc thetao_yearly_CNRM-CM6-1_historical_r${MEMBER}i1p1f2_gn_187501-189912.nc "${TMP_DIR}/thetao_yearly_CNRM-CM6-1_historical_r${MEMBER}i1p1f2_gn_185001-189912.nc"
fi

if [ -f thetao_yearly_CNRM-CM6-1_historical_r${MEMBER}i1p1f2_gn_190001-192412.nc ] && [ -f thetao_yearly_CNRM-CM6-1_historical_r${MEMBER}i1p1f2_gn_192501-194912.nc ]; then
    cdo mergetime thetao_yearly_CNRM-CM6-1_historical_r${MEMBER}i1p1f2_gn_190001-192412.nc thetao_yearly_CNRM-CM6-1_historical_r${MEMBER}i1p1f2_gn_192501-194912.nc "${TMP_DIR}/thetao_yearly_CNRM-CM6-1_historical_r${MEMBER}i1p1f2_gn_190001-194912.nc"
fi

if [ -f thetao_yearly_CNRM-CM6-1_historical_r${MEMBER}i1p1f2_gn_195001-197412.nc ] && [ -f thetao_yearly_CNRM-CM6-1_historical_r${MEMBER}i1p1f2_gn_197501-199912.nc ]; then
    cdo mergetime thetao_yearly_CNRM-CM6-1_historical_r${MEMBER}i1p1f2_gn_195001-197412.nc thetao_yearly_CNRM-CM6-1_historical_r${MEMBER}i1p1f2_gn_197501-199912.nc "${TMP_DIR}/thetao_yearly_CNRM-CM6-1_historical_r${MEMBER}i1p1f2_gn_195001-199912.nc"
fi

if [ -f thetao_yearly_CNRM-CM6-1_historical_r${MEMBER}i1p1f2_gn_200001-201412.nc ]; then
    cp thetao_yearly_CNRM-CM6-1_historical_r${MEMBER}i1p1f2_gn_200001-201412.nc "${TMP_DIR}/thetao_yearly_CNRM-CM6-1_historical_r${MEMBER}i1p1f2_gn_200001-201412.nc"
fi

# Remove original 25-year thetao files from BASE_DIR to avoid confusion
echo "Removing original 25-year thetao files"
rm -f thetao_yearly_CNRM-CM6-1_historical_r${MEMBER}i1p1f2_gn_185001-187412.nc thetao_yearly_CNRM-CM6-1_historical_r${MEMBER}i1p1f2_gn_187501-189912.nc \
      thetao_yearly_CNRM-CM6-1_historical_r${MEMBER}i1p1f2_gn_190001-192412.nc thetao_yearly_CNRM-CM6-1_historical_r${MEMBER}i1p1f2_gn_192501-194912.nc \
      thetao_yearly_CNRM-CM6-1_historical_r${MEMBER}i1p1f2_gn_195001-197412.nc thetao_yearly_CNRM-CM6-1_historical_r${MEMBER}i1p1f2_gn_197501-199912.nc

# --------------------------------------------------
# Step 2: convert merged thetao files from °C to K
# --------------------------------------------------

echo "Converting merged thetao files to Kelvin"

cd "${TMP_DIR}" || exit 1

for file in thetao_yearly_*.nc; do
    [ -e "$file" ] || continue

    out="${file/thetao/thetaoK}"

    if [ ! -f "$out" ]; then
        echo "Converting $file -> $out"
        cdo addc,273.15 "$file" "$out"
    fi
done

# --------------------------------------------------
# Step 3: compute 3D OHC and full-column 2D OHC
# --------------------------------------------------
for file in thetaoK_yearly_CNRM-CM6-1_historical_r${MEMBER}i1p1f2_gn_*.nc; do
    [ -e "$file" ] || continue

    period=$(basename "$file" | grep -oE '[0-9]{6}-[0-9]{6}') 
    thkcello_file="../thkcello_yearly_CNRM-CM6-1_historical_r${MEMBER}i1p1f2_gn_${period}.nc"

    if [ ! -f "$thkcello_file" ]; then
        echo "Missing thickness file for period ${period}: $thkcello_file"
        continue
    fi

    heatc3D_file="ohc_3D_yearly_${period}_r${MEMBER}.nc"
    heatc2D_file="ohc_2D_yearly_${period}_r${MEMBER}.nc"

    echo "--------------------------------------"
    echo "Processing period ${period}, member r${MEMBER}"
    echo "thetao file   : $file"
    echo "thkcello file : $thkcello_file"
    echo "--------------------------------------"

    # heatc3D = Cp * rho * T * dz
    # Units: J/m²
    cdo mulc,"$CPRHO" -mul "$file" "$thkcello_file" "$heatc3D_file"

    # full-column vertical integral
    # Units remain J/m²
    cdo vertsum "$heatc3D_file" "$heatc2D_file"
done

echo "Full-depth OHC completed for member r${MEMBER}"
