#!/bin/bash

set -euo pipefail

# ------------------------------------------------------------------
# Build historical AMOC members from CNRM-CM6-1 historical msftyz
#
# Output:
#   amoc_hist_r2.nc ... amoc_hist_r40.nc
# ------------------------------------------------------------------

EDIR_HIST="/cnrm/cmip/CMIP6/CMIP/CNRM-CERFACS/CNRM-CM6-1/historical"
OUT_DIR="/cnrm/ioga/Users/seguy/hist_ens/amoc_members"

mkdir -p "${OUT_DIR}"

python3 - <<PY
from pathlib import Path
import numpy as np
import xarray as xr
from dask.diagnostics import ProgressBar

edir_hist = Path(r"${EDIR_HIST}")
out_dir = Path(r"${OUT_DIR}")

var = "msftyz"
rho = 1026.0
echelle_amoc = 1e6

for member in range(2, 5):
    print(f"Processing historical AMOC r{member}")

    pattern = (
        edir_hist
        / f"r{member}i1p1f2"
        / "Omon"
        / "msftyz"
        / "gn"
        / "latest"
        / "*.nc"
    )

    fdata = xr.open_mfdataset(
        str(pattern),
        chunks={"time": 120},
        parallel=True,
        use_cftime=True,
    )

    vdata = fdata[var].isel(basin=1).groupby("time.year").mean(dim="time")

    if "j-mean" in vdata.dims:
        vdata = vdata.rename({"j-mean": "y"})

    # Historical should cover 1850-2014
    vdata = vdata.sel(year=slice(1850, 2014)).copy()

    if vdata.sizes["year"] != 165:
        raise ValueError(
            f"r{member}: expected 165 years, got {vdata.sizes['year']}"
        )

    # 214 = 45°N
    vmsf_45n = vdata.isel(y=214)
    vmsf_upper = vmsf_45n.sel(lev=slice(200, 2500))

    amoc = vmsf_upper.max("lev") / (rho * echelle_amoc)
    amoc = xr.where(amoc == 0, np.nan, amoc)

    ds_amoc = xr.Dataset({"amoc": amoc})

    output_path = out_dir / f"amoc_hist_r{member}.nc"
    print(f"Saving {output_path}")

    with ProgressBar():
        ds_amoc.to_netcdf(
            path=output_path,
            format="NETCDF3_CLASSIC",
        )

print("All historical AMOC members created.")
PY
