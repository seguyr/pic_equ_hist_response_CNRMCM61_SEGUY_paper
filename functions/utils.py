
"""
Utility functions for loading and processing climate diagnostics
used in the manuscript figures.
"""

from pathlib import Path
import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeature

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
N_BOOT = 1000
QT_INF = 0.05
QT_SUP = 0.95
seed = 0

COLORS = {
    "orange_dark": "#d95f02",
    "orange_light": "#fdd0a2",
    "teal_dark": "#1b9e77",
    "teal_light": "#a6dba0",
    "red_dark": "#b2182b",
    "red_light": "#f4a6b3",
    "blue_dark": "#2166ac",
    "blue_light": "#92c5de",
}

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

def load_integrated_ohc(dataset_type):
    if dataset_type == "hist_tot":
        filepath = DIR_HIST_TOT / "ohc_2D_hist_tot.nc"
    elif dataset_type == "hist_3000":
        filepath = DIR_HIST_3000 / "ohc_2D_hist_3000.nc"
    elif dataset_type == "pic_tot":
        filepath = DIR_PIC_TOT / "ohc_2D_pic_tot.nc"
    elif dataset_type == "pic_3000":
        filepath = DIR_PIC_3000 / "ohc_2D_pic_3000.nc"
    ds = xr.open_dataset(filepath)
    varname = list(ds.data_vars)[0]
    area_oce = xr.open_dataset(AREACELLO_FILE)["areacello"]
    return (ds[varname] * area_oce).sum(dim=("x", "y"))

def nemo_lon_lat(area_oce):
    """Return NEMO longitude/latitude with longitude discontinuity fixed."""
    area = area_oce.copy()
    after_discont = ~(area.lon.diff("x", label="upper") > 0).cumprod("x").astype(bool)
    area.lon[:, 1:] = area.lon[:, 1:] + after_discont * 360.0
    return area.lon, area.lat

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

def load_ohc_2d(dataset_type, layer="tot"):
    """
    Load one preprocessed 2D OHC ensemble field.
    """
    if dataset_type == "hist_tot":
        base_dir = DIR_HIST_TOT
        prefix = "ohc_2D_hist_tot"
    elif dataset_type == "hist_3000":
        base_dir = DIR_HIST_3000
        prefix = "ohc_2D_hist_3000"
    elif dataset_type == "pic_tot":
        base_dir = DIR_PIC_TOT
        prefix = "ohc_2D_pic_tot"
    elif dataset_type == "pic_3000":
        base_dir = DIR_PIC_3000
        prefix = "ohc_2D_pic_3000"
    else:
        raise ValueError(f"Unknown dataset_type: {dataset_type}")
    if layer == "tot":
        filepath = base_dir / f"{prefix}.nc"
    elif layer in ["0-300", "300-2000", "2000-bottom"]:
        filepath = base_dir / f"{prefix}_{layer}.nc"
    else:
        raise ValueError(f"Unknown layer: {layer}")
    ds = xr.open_dataset(filepath)
    return ds[list(ds.data_vars)[0]]

def time_matching(hist: xr.DataArray, pic: xr.DataArray) -> xr.DataArray:
    pic = pic.assign_coords(time=hist["time"])
    return hist - pic

def anomalies(arr: xr.DataArray, ref_slice=slice(0, 50)) -> xr.DataArray:
    ref = arr.isel(time=ref_slice).mean(dim="time")
    return arr - ref

def gain(arr: xr.DataArray, gain_slice=slice(145, 165)) -> xr.DataArray:
    return arr.isel(time=gain_slice).mean(dim="time")

def bootstrap(arr):
    """
    Bootstrap confidence intervals across the ensemble dimension.
    
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
        size=(N_BOOT, n_ens),
        replace=True,
    )
    sampled_da = xr.DataArray(sampled_ens, dims=("boot", "ensemble"))
    resampled = arr.sel(ensemble=sampled_da)
    boot_means = resampled.mean(dim="ensemble")
    boot_means = boot_means.chunk(dict(boot=-1))
    q_inf = boot_means.quantile(QT_INF, dim="boot")
    q_sup = boot_means.quantile(QT_SUP, dim="boot")
    mean = boot_means.mean(dim="boot")
    sigma = boot_means.std(dim="boot")
    p_left = (boot_means <= 0).mean(dim="boot")
    p_right = (boot_means >= 0).mean(dim="boot")
    p_value = 2 * np.minimum(p_left, p_right)
    return xr.concat(
        [q_inf, mean, q_sup, sigma, p_value],
        dim=xr.IndexVariable("stats", STATS_COORD),
    )

def bootstrap_2(arr1, arr2):
    """
    Bootstrap confidence intervals for the difference between two ensembles:
    mean(arr2) - mean(arr1)
    """
    if np.all(np.isnan(arr1)) or np.all(np.isnan(arr2)):
        out = xr.full_like(arr1.isel(ensemble=0, drop=True), np.nan).expand_dims(stats=STATS_COORD)
        return out.assign_coords(stats=STATS_COORD)
    rng = np.random.default_rng(seed)
    n_ens1 = arr1.sizes["ensemble"]
    ens_values1 = arr1.ensemble.values
    n_ens2 = arr2.sizes["ensemble"]
    ens_values2 = arr2.ensemble.values
    sampled_ens_1 = rng.choice(ens_values1, size=(N_BOOT, n_ens1), replace=True)
    sampled_ens_2 = rng.choice(ens_values2, size=(N_BOOT, n_ens2), replace=True)
    sampled_da_1 = xr.DataArray(sampled_ens_1, dims=("boot", "ensemble"))
    sampled_da_2 = xr.DataArray(sampled_ens_2, dims=("boot", "ensemble"))
    resampled_1 = arr1.sel(ensemble=sampled_da_1)
    resampled_2 = arr2.sel(ensemble=sampled_da_2)
    boot_means = resampled_2.mean(dim="ensemble") - resampled_1.mean(dim="ensemble")
    q_inf = boot_means.quantile(QT_INF, dim="boot")
    q_sup = boot_means.quantile(QT_SUP, dim="boot")
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

def boot(arr: xr.DataArray):
    """Run bootstrap and return a formatted xarray.DataArray."""
    result = bootstrap(arr)
    return to_bootstrap_xarray(result, arr)

def boot_diff(arr1: xr.DataArray, arr2: xr.DataArray):
    """Run bootstrap_2 and return a formatted xarray.DataArray."""
    result = bootstrap_2(arr1, arr2)
    return to_bootstrap_xarray(result, arr1)

def get_stats(arr: xr.DataArray):
    lower = arr.sel(stats="lower")
    mean = arr.sel(stats="mean")
    upper = arr.sel(stats="upper")
    return lower, mean, upper

def get_ci(arr: xr.DataArray):
    lower, mean, upper = get_stats(arr)
    return mean, lower, upper

def get_scalar_stats(arr: xr.DataArray):
    lower, mean, upper = get_stats(arr)
    return float(lower.values), float(mean.values), float(upper.values)

def remove_map_outline(ax):
    """Enlève complètement le contour noir (outline) de la projection."""
    if hasattr(ax, "outline_patch"):
        ax.outline_patch.set_visible(False)
    if "geo" in ax.spines:
        ax.spines["geo"].set_visible(False)
    ax.patch.set_edgecolor("none")
    ax.patch.set_linewidth(0)

def plot_panel(ax, ds, title, label, cmap, norm):
    """
    Plot one spatial panel from a bootstrap DataArray.
    """
    low, mean, up = get_stats(ds)
    area_oce = load_area_ocean()
    lon, lat = nemo_lon_lat(area_oce)
    ax.set_global()
    ax.coastlines(linewidth=0.55)
    ax.add_feature(cfeature.BORDERS, linewidth=0.35, linestyle=":")
    cf = ax.pcolormesh(
        lon,
        lat,
        mean,
        cmap=cmap,
        norm=norm,
        transform=ccrs.PlateCarree(),
        shading="auto",
    )
    # Hatch where confidence interval includes zero
    mask_non_sig = np.ma.masked_where((low >= 0) | (up <= 0), mean)
    hatch = ax.contourf(
        lon,
        lat,
        mask_non_sig,
        hatches=["///"],
        colors="none",
        transform=ccrs.PlateCarree(),
        zorder=2,
    )
    for coll in hatch.collections:
        coll.set_edgecolor("lightgray")
        coll.set_linewidth(0.0)
    ax.set_title(title, fontsize=16, pad=8, fontweight="bold")
    ax.text(
        0.02,
        0.96,
        label,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=18,
        fontweight="bold",
    )
    remove_map_outline(ax)
    return cf

