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

MEMBER="${1:?Please provide a member number, e.g. 1 or 25}"

OUT_DIR="/cnrm/ioga/Users/seguy/hist_ens/hist_${MEMBER}"
THKCELLO_DIR="/cnrm/cmip6/CMIP/CNRM-CERFACS/CNRM-CM6-1/historical/r${MEMBER}i1p1f2/Omon/thkcello/gn/latest"
THETAO_DIR="/cnrm/cmip6/CMIP/CNRM-CERFACS/CNRM-CM6-1/historical/r${MEMBER}i1p1f2/Omon/thetao/gn/latest"

mkdir -p "${OUT_DIR}"

echo "=== Processing member r${MEMBER}i1p1f2 ==="
echo "Output directory: ${OUT_DIR}"

# --- thkcello ---
cd "${THKCELLO_DIR}" || exit 1
for filet in thkcello_*.nc; do
    [ -e "$filet" ] || continue
    outfile="${OUT_DIR}/${filet/Omon/yearly}"
    echo "Processing thkcello file: $filet"
    cdo -L yearmean "$filet" "$outfile"
done

# --- thetao ---
cd "${THETAO_DIR}" || exit 1
for filet in thetao_*.nc; do
    [ -e "$filet" ] || continue
    outfile="${OUT_DIR}/${filet/Omon/yearly}"
    echo "Processing thetao file: $filet"
    cdo -L yearmean "$filet" "$outfile"
done

echo "Yearly means completed for member ${MEMBER}"
