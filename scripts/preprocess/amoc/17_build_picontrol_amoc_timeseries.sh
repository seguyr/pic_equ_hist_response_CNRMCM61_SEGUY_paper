#!/bin/bash


set -euo pipefail

# ------------------------------------------------------------------
# Build full piControl AMOC time series from msftyz
#
# Output:
#   amoc_picontrol.nc
#
# AMOC definition:
#   maximum Atlantic meridional overturning streamfunction
#   at 45°N, between 200 and 2500 m depth
#
# Units:
#   Sverdrup (Sv)
# ------------------------------------------------------------------

MSFTYZ_GLOB="/cnrm/cmip/CMIP6/CMIP/CNRM-CERFACS/CNRM-CM6-1/piControl/r1i1p1f2/Omon/msftyz/gn/latest/*.nc"
OUT_DIR="/cnrm/ioga/Users/seguy/pic_ens/pic"
OUT_FILE="${OUT_DIR}/amoc_picontrol.nc"

mkdir -p "${OUT_DIR}"

python3 - <<PY
from pathlib import Path
import numpy as np
import xarray as xr
from dask.diagnostics import ProgressBar

msftyz_glob = r"${MSFTYZ_GLOB}"
out_file = Path(r"${OUT_FILE}")

var = "msftyz"
rho = 1026.0
echelle_amoc = 1e6

print("Opening piControl msftyz...")
fdata = xr.open_mfdataset(
    msftyz_glob,
    chunks={"time": 120},
    parallel=True,
    use_cftime=True,
)

print("Computing yearly mean streamfunction...")
vdata = fdata[var].isel(basin=1).groupby("time.year").mean(dim="time")

if "j-mean" in vdata.dims:
    vdata = vdata.rename({"j-mean": "y"})

print("Calcul de l'AMOC : maximum de la fonction de courant méridienne en Atlantique, à 45°N, entre 200 et 2500 m")

# 214 = 45°N in your grid convention
vmsf_45n = vdata.isel(y=214)

# Upper cell between 200 and 2500 m
vmsf_upper = vmsf_45n.sel(lev=slice(200, 2500))

# Maximum over depth, converted to Sv
amoc = vmsf_upper.max("lev") / (rho * echelle_amoc)

# Replace zeros by NaN if needed
amoc = xr.where(amoc == 0, np.nan, amoc)

ds_amoc = xr.Dataset({"amoc": amoc})

print(f"Saving {out_file}")
with ProgressBar():
    ds_amoc.to_netcdf(
        path=out_file,
        format="NETCDF3_CLASSIC",
    )

print("Done.")
PY
