
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

STATS_COORD = ["lower", "mean", "upper", "sigma", "p_value"]

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
    area_oce = xr.open_dataset(AREACELLO_FILE)["areacello"]
    return (ds[varname] * area_oce).sum(dim=("x", "y"))


def load_branching_years():
    """Return branching years for both ensembles."""
    hist_starty = xr.open_dataset(HIST_START_YEARS_FILE)["year"]
    starty_tot = hist_starty[0:29].values.tolist()
    starty_3000 = hist_starty[29:39].values.tolist()
    return starty_tot, starty_3000


def _first_data_var(ds: xr.Dataset) -> xr.DataArray:
    """Return the first data variable of a Dataset."""
    return ds[list(ds.data_vars)[0]]


def load_ohc_2d_ensembles():
    hist_tot = _first_data_var(xr.open_dataset(DIR_HIST_TOT / "ohc_2D_hist_tot.nc"))
    hist_3000 = _first_data_var(xr.open_dataset(DIR_HIST_3000 / "ohc_2D_hist_3000.nc"))
    pic_tot = _first_data_var(xr.open_dataset(DIR_PIC_TOT / "ohc_2D_pic_tot.nc"))
    pic_3000 = _first_data_var(xr.open_dataset(DIR_PIC_3000 / "ohc_2D_pic_3000.nc"))
    return hist_tot, hist_3000, pic_tot, pic_3000


def time_matching(hist: xr.DataArray, pic: xr.DataArray) -> xr.DataArray:
    pic = pic.assign_coords(time=hist["time"])
    return hist - pic


def anomalies(arr: xr.DataArray, ref_slice=slice(0, 50)) -> xr.DataArray:
    ref = arr.isel(time=ref_slice).mean(dim="time")
    return arr - ref


def gain(arr: xr.DataArray, gain_slice=slice(145, 165)) -> xr.DataArray:
    return arr.isel(time=gain_slice).mean(dim="time")


# ============================================================
# Bootstrap functions
# ============================================================

def bootstrap(arr, n_boot=1000, qt_inf=0.05, qt_sup=0.95, seed=0):
    """
    Bootstrap confidence intervals across the ensemble dimension.

    Parameters
    ----------
    arr : xr.DataArray
        Input array with an 'ensemble' dimension.
        Can also contain additional dimensions such as time, x, y.
    n_boot : int
        Number of bootstrap resamples.
    qt_inf, qt_sup : float
        Lower and upper quantiles.
    seed : int
        Random seed.

    Returns
    -------
    xr.DataArray
        DataArray with dimensions ('stats', ...)
    """
    if np.all(np.isnan(arr)):
        out_shape = tuple(arr.sizes[d] for d in arr.dims if d != "ensemble")
        out_dims = tuple(d for d in arr.dims if d != "ensemble")
        out_coords = {d: arr.coords[d] for d in out_dims}
        out = xr.full_like(arr.isel(ensemble=0, drop=True), np.nan).expand_dims(stats=STATS_COORD)
        return out.assign_coords(stats=STATS_COORD)
    rng = np.random.default_rng(seed)
    n_ens = arr.sizes["ensemble"]
    ens_values = arr.ensemble.values
    sampled_ens = rng.choice(
        ens_values,
        size=(n_boot, n_ens),
        replace=True,
    )
    sampled_da = xr.DataArray(sampled_ens, dims=("boot", "ensemble"))
    resampled = arr.sel(ensemble=sampled_da)
    boot_means = resampled.mean(dim="ensemble")
    boot_means = boot_means.chunk(dict(boot=-1))
    q_inf = boot_means.quantile(qt_inf, dim="boot")
    q_sup = boot_means.quantile(qt_sup, dim="boot")
    mean = boot_means.mean(dim="boot")
    sigma = boot_means.std(dim="boot")
    p_left = (boot_means <= 0).mean(dim="boot")
    p_right = (boot_means >= 0).mean(dim="boot")
    p_value = 2 * np.minimum(p_left, p_right)
    return xr.concat(
        [q_inf, mean, q_sup, sigma, p_value],
        dim=xr.IndexVariable("stats", STATS_COORD),
    )


def bootstrap_2(arr1, arr2, n_boot=1000, qt_inf=0.025, qt_sup=0.975, seed=0):
    """
    Bootstrap confidence intervals for the difference between two ensembles:
    mean(arr2) - mean(arr1)

    Parameters
    ----------
    arr1, arr2 : xr.DataArray
        Input arrays with an 'ensemble' dimension.
    n_boot : int
        Number of bootstrap resamples.
    qt_inf, qt_sup : float
        Lower and upper quantiles.
    seed : int
        Random seed.

    Returns
    -------
    xr.DataArray
        DataArray with dimensions ('stats', ...)
    """
    if np.all(np.isnan(arr1)) or np.all(np.isnan(arr2)):
        out = xr.full_like(arr1.isel(ensemble=0, drop=True), np.nan).expand_dims(stats=STATS_COORD)
        return out.assign_coords(stats=STATS_COORD)
    rng = np.random.default_rng(seed)
    n_ens1 = arr1.sizes["ensemble"]
    ens_values1 = arr1.ensemble.values
    n_ens2 = arr2.sizes["ensemble"]
    ens_values2 = arr2.ensemble.values
    sampled_ens_1 = rng.choice(ens_values1, size=(n_boot, n_ens1), replace=True)
    sampled_ens_2 = rng.choice(ens_values2, size=(n_boot, n_ens2), replace=True)
    sampled_da_1 = xr.DataArray(sampled_ens_1, dims=("boot", "ensemble"))
    sampled_da_2 = xr.DataArray(sampled_ens_2, dims=("boot", "ensemble"))
    resampled_1 = arr1.sel(ensemble=sampled_da_1)
    resampled_2 = arr2.sel(ensemble=sampled_da_2)
    boot_means = resampled_2.mean(dim="ensemble") - resampled_1.mean(dim="ensemble")
    q_inf = boot_means.quantile(qt_inf, dim="boot")
    q_sup = boot_means.quantile(qt_sup, dim="boot")
    mean = boot_means.mean(dim="boot")
    sigma = boot_means.std(dim="boot")
    p_left = (boot_means <= 0).mean(dim="boot")
    p_right = (boot_means >= 0).mean(dim="boot")
    p_value = 2 * np.minimum(p_left, p_right)
    return xr.concat(
        [q_inf, mean, q_sup, sigma, p_value],
        dim=xr.IndexVariable("stats", STATS_COORD),
    )

def to_bootstrap_xarray(result, template: xr.DataArray) -> xr.DataArray:
    out_dims = tuple(d for d in template.dims if d != "ensemble")
    out_coords = {d: template.coords[d] for d in out_dims}
    return xr.DataArray(
        result,
        dims=("stats",) + out_dims,
        coords={"stats": STATS_COORD, **out_coords},
    )

def boot(arr: xr.DataArray, n_boot=1000, qt_inf=0.025, qt_sup=0.975, seed=0):
    """Run bootstrap and return a formatted xarray.DataArray."""
    result = bootstrap(arr, n_boot=n_boot, qt_inf=qt_inf, qt_sup=qt_sup, seed=seed)
    return to_bootstrap_xarray(result, arr)


def boot_diff(arr1: xr.DataArray, arr2: xr.DataArray, n_boot=1000, qt_inf=0.025, qt_sup=0.975, seed=0):
    """Run bootstrap_2 and return a formatted xarray.DataArray."""
    result = bootstrap_2(arr1, arr2, n_boot=n_boot, qt_inf=qt_inf, qt_sup=qt_sup, seed=seed)
    return to_bootstrap_xarray(result, arr1)
