#!/bin/bash
#
#set -euo pipefail

# ------------------------------------------------------------------
# Build 39 piControl pseudo-members of AMOC from msftyz
#
# This script:
#   1. reads branching years from data/metadata/hist_start_years.nc
#   2. extracts 165-year piControl segments for each pseudo-member
#   3. reindexes time to 1850-2014
#   4. computes AMOC at 45°N, between 200 and 2500 m
#   5. converts the result to Sverdrup (Sv)
#
# Output files:
#   amoc_pic_r1.nc ... amoc_pic_r39.nc
# ------------------------------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(dirname "$(dirname "$SCRIPT_DIR")")"

BRANCH_FILE="${REPO_ROOT}/data/metadata/hist_start_years.nc"
OUT_DIR="/cnrm/ioga/Users/seguyr/pic_ens/pic_amoc"
MSFTYZ_GLOB="/cnrm/cmip/CMIP6/CMIP/CNRM-CERFACS/CNRM-CM6-1/piControl/r1i1p1f2/Omon/msftyz/gn/latest/*.nc"

mkdir -p "${OUT_DIR}"

if [ ! -f "${BRANCH_FILE}" ]; then
    echo "Missing metadata file: ${BRANCH_FILE}"
    exit 1
fi

python3 - <<PY
from pathlib import Path
import numpy as np
import xarray as xr
from dask.diagnostics import ProgressBar

branch_file = Path(r"${BRANCH_FILE}")
out_dir = Path(r"${OUT_DIR}")
msftyz_glob = r"${MSFTYZ_GLOB}"

var = "msftyz"
rho = 1026.0
echelle_amoc = 1e6

print("Opening piControl msftyz...")
fdata = xr.open_mfdataset(
    msftyz_glob,
    parallel=True,
    use_cftime=True,
)

print("Computing yearly mean streamfunction...")
vdata = fdata[var].isel(basin=1).groupby("time.year").mean(dim="time")

if "j-mean" in vdata.dims:
    vdata = vdata.rename({"j-mean": "y"})

start_years = xr.open_dataset(branch_file)["year"].values.tolist()
start_years = [int(y) for y in start_years]

print(f"Number of pseudo-members: {len(start_years)}")

for i, iniyear in enumerate(start_years, start=1):
    endyear = iniyear + 164
    print(f"Processing pic_r{i}: {iniyear}-{endyear}")

    vdata_1d = vdata.sel(year=slice(iniyear, endyear)).copy()

    if vdata_1d.sizes["year"] != 165:
        raise ValueError(
            f"pic_r{i}: expected 165 years, got {vdata_1d.sizes['year']} "
            f"for slice {iniyear}-{endyear}"
        )

    vdata_1d = vdata_1d.assign_coords(year=np.arange(1850, 2014))

    # 214 = 45°N in your grid convention
    vmsf_45n = vdata_1d.isel(y=214)

    # Upper cell between 200 and 2500 m
    vmsf_upper = vmsf_45n.sel(lev=slice(200, 2500))

    # Maximum over depth, then convert to Sv
    amoc = vmsf_upper.max("lev") / (rho * echelle_amoc)

    amoc = xr.where(amoc == 0, np.nan, amoc)

    ds_amoc = xr.Dataset({"amoc": amoc})

    output_path = out_dir / f"amoc_pic_r{i}.nc"
    print(f"Saving {output_path}")

    with ProgressBar():
        ds_amoc.to_netcdf(
            path=output_path,
            format="NETCDF3_CLASSIC",
        )

print("All piControl AMOC pseudo-members created.")
PY
