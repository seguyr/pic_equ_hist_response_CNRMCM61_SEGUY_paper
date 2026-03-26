
"""
Utility functions for loading and processing climate diagnostics
used in the manuscript figures.
"""


from pathlib import Path
import numpy as np
import xarray as xr


def load_area_ocean(area_dir: Path) -> xr.DataArray:
    """Load ocean grid-cell area."""
    ds = xr.open_mfdataset(str(area_dir / "*.nc"))
    return ds["areacello"]


def load_pic_ohc(pic_dir, cp, rho, area_tot, scale):
    """Load piControl global OHC."""
    ds = xr.open_mfdataset(str(pic_dir / "heatc1D*.nc"))

    ohc = (
        ds.thetao.isel(lon=0, lat=0)
        * cp * rho * area_tot
        / scale
    )

    return ohc.assign_coords(time=np.arange(3000))


def load_pic_amoc(amoc_dir, rho, scale):
    """Load piControl AMOC."""
    ds = xr.open_mfdataset(str(amoc_dir / "amoc_pic_gb_3000y*.nc"))

    amoc = ds.msftyz / (rho * scale)

    return amoc.assign_coords(year=np.arange(3001))


def load_integrated_ohc(filepath, area_oce):
    """Horizontally integrate OHC field."""
    ds = xr.open_mfdataset(str(filepath))

    return (ds * area_oce).sum(dim=("x", "y")).__xarray_dataarray_variable__


def load_branching_years(filepath):
    """Return branching years for both ensembles."""
    hist_starty = xr.open_dataset(filepath)["year"]

    starty_tot = hist_starty[0:29].values.tolist()
    starty_3000 = hist_starty[29:39].values.tolist()

    return starty_tot, starty_3000
