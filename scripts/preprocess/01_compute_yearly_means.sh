#!/bin/bash

# ------------------------------------------------------------------
# Compute yearly means from CMIP6 monthly ocean variables
#
# This script converts monthly CMIP6 ocean temperature (thetao) and
# layer thickness (thkcello) fields from the CNRM-CM6-1 historical
# simulations into yearly means using CDO.
#
# Input files are processed separately for each ensemble member and
# written as yearly NetCDF time series covering 1850–2014.
#
# These yearly fields are used as input for the ocean heat content
# (OHC) computation performed in subsequent preprocessing steps.
# ------------------------------------------------------------------


set -euo pipefail

# Compute yearly means for one historical ensemble member
# Usage:
#   bash scripts/preprocess/01_compute_yearly_means.sh 1
#   bash scripts/preprocess/01_compute_yearly_means.sh 25

MEMBER="${1:?Please provide a member number, e.g. 1 or 25}"

REPO="/cnrm/ioga/Users/seguy/hist_ens/hist_${MEMBER}"
THKCELLO_DIR="/cnrm/cmip6/CMIP/CNRM-CERFACS/CNRM-CM6-1/historical/r${MEMBER}i1p1f2/Omon/thkcello/gn/latest"
THETAO_DIR="/cnrm/cmip6/CMIP/CNRM-CERFACS/CNRM-CM6-1/historical/r${MEMBER}i1p1f2/Omon/thetao/gn/latest"

mkdir -p "${REPO}"

echo "=== Processing member r${MEMBER}i1p1f2 ==="
echo "Output directory: ${REPO}"

# --- thkcello ---
cd "${THKCELLO_DIR}" || exit 1
for filet in thkcello_*.nc; do
    [ -e "$filet" ] || continue
    file="${filet//@/}"
    echo "Processing thkcello file: $filet"
    ncks -3 "$filet" "${REPO}/${file}"
    cdo yearavg "${REPO}/${file}" "${REPO}/${file/Omon/yearly}"
    rm -f "${REPO}/${file}"
done

# --- thetao ---
cd "${THETAO_DIR}" || exit 1
for filet in thetao_*.nc; do
    [ -e "$filet" ] || continue
    file="${filet//@/}"
    echo "Processing thetao file: $filet"
    ncks -3 "$filet" "${REPO}/${file}"
    cdo yearavg "${REPO}/${file}" "${REPO}/${file/Omon/yearly}"
    rm -f "${REPO}/${file}"
done

echo "Yearly means completed for member ${MEMBER}"
