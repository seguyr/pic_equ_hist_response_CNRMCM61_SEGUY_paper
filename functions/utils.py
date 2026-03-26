
"""
Utility functions for loading and processing climate diagnostics
used in the manuscript figures.
"""

from pathlib import Path
import numpy as np
import xarray as xr

REPO_ROOT = Path(__file__).resolve().parents[1]

METADATA_DIR = REPO_ROOT / "data" / "metadata"

HIST_START_YEARS_FILE = METADATA_DIR / "hist_start_years.nc"
AREACELLO_FILE = METADATA_DIR / "areacello.nc"
MASK_FILE = METADATA_DIR / "mask.nc"

DATA_ROOT = Path("/cnrm/ioga/Users/seguyr")

DIR_HIST_TOT = DATA_ROOT / "hist_ens/hist_tot"
DIR_HIST_3000 = DATA_ROOT / "hist_ens/hist_3000"

DIR_PIC_TOT = DATA_ROOT / "pic_ens/pic_tot"
DIR_PIC_3000 = DATA_ROOT / "pic_ens/pic_3000"

DIR_PIC = DATA_ROOT / "pic_ens/pic"

DIR_AMOC = Path("/cnrm/ioga/Users/seguy/hist_ensbl/pic_global")


def load_area_ocean():
    """Load ocean grid-cell area from repository metadata."""
    ds = xr.open_dataset(AREACELLO_FILE)
    return ds["areacello"]


def load_pic_ohc(pic_dir, scale):
    """Load piControl global OHC."""
    ds = xr.open_mfdataset(str(pic_dir / "heatc1D*.nc"), use_cftime=True)

    ohc = (
        ds.thetao / scale
    )

    return ohc.assign_coords(time=np.arange(3000))


def load_pic_amoc(amoc_dir, rho, scale):
    """Load piControl AMOC."""
    ds = xr.open_mfdataset(str(amoc_dir / "amoc_pic_gb_3000y*.nc"), use_cftime=True)
    
    amoc = ds.amoc

    return amoc.assign_coords(year=np.arange(3001))


def load_integrated_ohc(filepath, area_oce):
    """Horizontally integrate OHC field."""
    ds = xr.open_mfdataset(str(filepath))

    return (ds * area_oce).sum(dim=("x", "y")).__xarray_dataarray_variable__


def load_branching_years():
    """Return branching years for both ensembles."""

    hist_starty = xr.open_dataset(HIST_START_YEARS_FILE)["year"]

    starty_tot = hist_starty[0:29].values.tolist()
    starty_3000 = hist_starty[29:39].values.tolist()

    return starty_tot, starty_3000
