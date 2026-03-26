
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


def load_area_ocean():
    """Load ocean grid-cell area from repository metadata."""
    ds = xr.open_dataset(AREACELLO_FILE)
    return ds["areacello"]


def load_pic_ohc(scale):
    """Load full piControl global OHC time series (scaled to ZJ)."""
    ds = xr.open_dataset(DIR_PIC / "ohc_1D_picontrol.nc")
    varname = list(ds.data_vars)[0]
    ohc = ds[varname] / scale
    return ohc.assign_coords(time=np.arange(ohc.sizes["time"]))


def load_pic_amoc():
    """Load full piControl AMOC time series (Sv)."""
    ds = xr.open_dataset(DIR_PIC / "amoc_picontrol.nc")
    amoc = ds["amoc"]
    return amoc.assign_coords(year=np.arange(amoc.sizes["year"]))


def load_integrated_ohc(dataset_type, area_oce):
    if dataset_type == "hist_tot":
        filepath = DIR_HIST_TOT / "ohc_2D_hist_tot.nc"
    elif dataset_type == "hist_3000":
        filepath = DIR_HIST_3000 / "ohc_2D_hist_3000.nc"
    ds = xr.open_dataset(filepath)
    varname = list(ds.data_vars)[0]
    return (ds[varname] * area_oce).sum(dim=("x", "y"))


def load_branching_years():
    """Return branching years for both ensembles."""
    hist_starty = xr.open_dataset(HIST_START_YEARS_FILE)["year"]
    starty_tot = hist_starty[0:29].values.tolist()
    starty_3000 = hist_starty[29:39].values.tolist()
    return starty_tot, starty_3000
