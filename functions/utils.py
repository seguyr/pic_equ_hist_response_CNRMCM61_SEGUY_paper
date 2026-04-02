
"""
Utility functions for loading and processing climate diagnostics
used in the manuscript figures.
"""

from pathlib import Path
import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import re
import dask
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors


REPO_ROOT = Path(__file__).resolve().parents[1]

METADATA_DIR = REPO_ROOT / "data" / "metadata"
PIC_GLOBAL_DIR = REPO_ROOT / "data" / "pic_global"
TAS_PIC_FILE = PIC_GLOBAL_DIR / "tas_pic_3000y.nc"
TOS_PIC_FILE = PIC_GLOBAL_DIR / "tos_pic_3000y.nc"
MOC_PIC_FILE = PIC_GLOBAL_DIR / "moc_pic_3000y.nc"
DOT_PIC_FILE = PIC_GLOBAL_DIR / "dot_pic_3000y.nc"
TOA_PIC_FILE = PIC_GLOBAL_DIR / "toa_pic_3000y.nc"
OHC_PIC_FILE = PIC_GLOBAL_DIR / "heatc1D_yearly_CNRM-CM6-1_piControl_r1i1p1f2_gn_185001-484912.nc"
HIST_START_YEARS_FILE = METADATA_DIR / "hist_start_years.nc"
AREACELLO_FILE = METADATA_DIR / "areacello_Ofx_CNRM-CM6-1_piControl_r1i1p1f2_gn.nc"
AREACELLA_FILE = METADATA_DIR / "areacella_fx_CNRM-CM6-1_piControl_r1i1p1f2_gr.nc"
MASK_FILE = METADATA_DIR / "gsubpolmsk.nc"


DIR_HIST_TOT = REPO_ROOT / "data" / "hist_tot"
DIR_HIST_3000 = REPO_ROOT / "data" / "hist_3000"
DIR_PIC_TOT = REPO_ROOT / "data" / "pic_tot"
DIR_PIC_3000 = REPO_ROOT / "data" / "pic_3000"

STATS_COORD = ["lower", "mean", "upper", "sigma", "p_value"]
N_BOOT = 1000
QT_INF = 0.05
QT_SUP = 0.95
seed = 0
CP =  3991.867 # J/(kg*K)
RHO = 1026 # kg/m³ 
AREA_TOT = 3.626975*10**14
ECHELLE_AMOC = 1e6


COLORS = {
    "orange_dark": "#e07b00",
    "orange_light": "#fdd87a",
    "teal_dark": "#005f73",
    "teal_light": "#9ae3d2",
    "red_dark": "#b2182b",
    "red_light": "#f4a582",
    "blue_dark": "#1b365d",
    "blue_light": "#a3bce2",
}
WHITE_FRAC = 0.01


def load_area_ocean():
    """Load ocean grid-cell area from repository metadata."""
    ds = xr.open_dataset(AREACELLO_FILE)
    return ds["areacello"]


def load_pic_ohc_1D():
    """Load full piControl global OHC time series."""
    ds = xr.open_dataset(OHC_PIC_FILE)
    ohc = ds.thetao.isel(lon=0).isel(lat=0)*CP*RHO*AREA_TOT
    return ohc.assign_coords(time=np.arange(ohc.sizes["time"]))

    
def load_2D_ohc(dataset_type, layer):
    if dataset_type == "hist_tot":
        filepath = DIR_HIST_TOT
    elif dataset_type == "hist_3000":
        filepath = DIR_HIST_3000
    elif dataset_type == "pic_tot":
        filepath = DIR_PIC_TOT
    elif dataset_type == "pic_3000":
        filepath = DIR_PIC_3000
    else:
        raise ValueError(f"Unknown dataset_type: {dataset_type}")
    files = sorted(Path(filepath).glob(f"*{layer}*.nc"), key=natural_key)
    ds = xr.open_mfdataset(
        files,
        combine="nested",
        concat_dim="ensemble"
    )
    return ds.__xarray_dataarray_variable__.assign_coords(time=np.arange(ds.sizes["time"]))


def load_pic_amoc():
    """Load full piControl AMOC time series (Sv)."""
    ds = xr.open_dataset(PIC_GLOBAL_DIR / "amoc_pic_gb_3000y.nc")
    amoc = ds.msftyz/(RHO*ECHELLE_AMOC)
    return amoc.assign_coords(year=np.arange(amoc.sizes["year"]))

def natural_key(path):
    return [int(text) if text.isdigit() else text.lower()
            for text in re.split(r'(\d+)', str(path))]

def load_integrated_ohc(dataset_type, layer):
    if dataset_type == "hist_tot":
        filepath = DIR_HIST_TOT
    elif dataset_type == "hist_3000":
        filepath = DIR_HIST_3000
    elif dataset_type == "pic_tot":
        filepath = DIR_PIC_TOT
    elif dataset_type == "pic_3000":
        filepath = DIR_PIC_3000
    else:
        raise ValueError(f"Unknown dataset_type: {dataset_type}")
    files = sorted(Path(filepath).glob(f"*{layer}*.nc"), key=natural_key)
    ds = xr.open_mfdataset(
        files,
        combine="nested",
        concat_dim="ensemble"
    )
    area_oce = xr.open_dataset(AREACELLO_FILE)["areacello"]
    return (ds.__xarray_dataarray_variable__ * area_oce).sum(dim=("x", "y"))

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
    

def load_ohc_2d(dataset_type, layer="all"):
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
    if layer == "all":
        filepath = base_dir / f"{prefix}.nc"
    elif layer in ["0-300", "300-2000", "2000-bottom"]:
        filepath = base_dir / f"{prefix}_{layer}.nc"
    else:
        raise ValueError(f"Unknown layer: {layer}")
    ds = xr.open_dataset(filepath)
    return ds[list(ds.data_vars)[0]]

def make_cmap_norm(vmin, vmax, white_frac=0.01, n=256):
    white_width = white_frac * (vmax - vmin)

    bounds = np.linspace(vmin, vmax, n)
    colors = plt.cm.RdBu_r(np.linspace(0, 1, n))

    i0_low = np.argmin(np.abs(bounds + white_width))
    i0_high = np.argmin(np.abs(bounds - white_width))
    colors[i0_low:i0_high] = [1, 1, 1, 1]

    cmap = mcolors.ListedColormap(colors)
    norm = mcolors.TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)
    return cmap, norm



def time_matching(hist: xr.DataArray, pic: xr.DataArray) -> xr.DataArray:
    pic = pic.assign_coords(time=hist["time"])
    return hist - pic

def anomalies(arr: xr.DataArray, ref_slice=slice(0, 50)) -> xr.DataArray:
    ref = arr.isel(time=ref_slice).mean(dim="time")
    return arr - ref

def gain(arr: xr.DataArray, gain_slice=slice(145, 165)) -> xr.DataArray:
    return arr.isel(time=gain_slice).mean(dim="time")

def bootstrap(arr, n_boot=N_BOOT, seed=0):
    if np.all(np.isnan(arr)):
        return np.array([np.nan, np.nan, np.nan, np.nan, np.nan], dtype=np.float32)

    # Générateur pseudo-aléatoire
    rng = np.random.default_rng(seed)

    n_ens = arr.sizes["ensemble"]
    ens_values = arr.ensemble.values

    print("sample")
    sampled_ens = rng.choice(
        ens_values,
        size=(n_boot, n_ens),
        replace=True
    )  # (n_boot, n_ens)

    sampled_da = xr.DataArray(sampled_ens, dims=("boot", "ensemble"))

    print("resample")
    resampled = arr.sel(ensemble=sampled_da)

    print("mean")
    boot_means = resampled.mean(dim="ensemble")
    boot_means = boot_means.chunk(dict(boot=-1))
    
    print("stats")
    q_inf = boot_means.quantile(QT_INF, dim="boot")
    q_sup = boot_means.quantile(QT_SUP, dim="boot")
    mean  = boot_means.mean(dim="boot")
    sigma = boot_means.std(dim="boot")

    p_left  = (boot_means <= 0).mean(dim="boot")
    p_right = (boot_means >= 0).mean(dim="boot")
    p_value = 2 * np.minimum(p_left, p_right)

    return np.array([q_inf, mean, q_sup, sigma, p_value], dtype=np.float32)


def bootstrap_2(arr1, arr2, n_boot = N_BOOT, seed=0):
    if np.all(np.isnan(arr1)) or np.all(np.isnan(arr2)):
        return np.array([np.nan, np.nan, np.nan, np.nan, np.nan], dtype=np.float32)

    rng = np.random.default_rng(seed)

    n_ens1 = arr1.sizes["ensemble"]
    ens_values1 = arr1.ensemble.values
    
    n_ens2 = arr2.sizes["ensemble"]
    ens_values2 = arr2.ensemble.values
    
    sampled_ens_1 = rng.choice(ens_values1, size=(n_boot, n_ens1), replace=True)
    sampled_ens_2 = rng.choice(ens_values2, size=(n_boot, n_ens2), replace=True)  # indépendant du 1

    sampled_da_1 = xr.DataArray(sampled_ens_1, dims=("boot", "ensemble"))
    sampled_da_2 = xr.DataArray(sampled_ens_2, dims=("boot", "ensemble"))

    resampled_1 = arr1.sel(ensemble=sampled_da_1)
    resampled_2 = arr2.sel(ensemble=sampled_da_2)

    boot_means = resampled_2.mean(dim="ensemble") - resampled_1.mean(dim="ensemble")

    q_inf = boot_means.quantile(QT_INF, dim="boot")
    q_sup = boot_means.quantile(QT_SUP, dim="boot")
    mean  = boot_means.mean(dim="boot")
    sigma = boot_means.std(dim="boot")

    p_left  = (boot_means <= 0).mean(dim="boot")
    p_right = (boot_means >= 0).mean(dim="boot")
    p_value = 2 * np.minimum(p_left, p_right)

    return np.array([q_inf, mean, q_sup, sigma, p_value], dtype=np.float32)


def get_stats(arr):
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
    
def hac_to_stats_da(low, mean, up):
    return xr.concat(
        [low,mean,up,
            xr.full_like(mean, np.nan),
            xr.full_like(mean, np.nan),
        ],
        dim=xr.IndexVariable("stats", STATS_COORD),
    )
    
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

def plot_panel_2(ax, title, panel_label, low, mean, up, cmap, norm_cmap):
    ax.set_global()
    area_oce = load_area_ocean()
    lon2d, lat2d = nemo_lon_lat(area_oce)
    ax.coastlines(linewidth=0.55)
    ax.add_feature(cfeature.BORDERS, linewidth=0.35, linestyle=':')

    m = ax.pcolormesh(
        lon2d, lat2d, mean,
        cmap=cmap, norm=norm_cmap,
        transform=ccrs.PlateCarree()
    )

    mask_non_sig = np.ma.masked_where((low >= 0) | (up <= 0), mean)
    hatch = ax.contourf(
        lon2d, lat2d, mask_non_sig,
        hatches=['///'],
        colors='none',
        transform=ccrs.PlateCarree(),
        zorder=2
    )

    
    # Style des hachures, compatible selon versions
    if hasattr(hatch, "collections"):
        for coll in hatch.collections:
            coll.set_edgecolor("gray")
            coll.set_linewidth(0.0)
    else:
        hatch.set_edgecolor("gray")
        hatch.set_linewidth(0.0)

    #ax.set_title(title, fontsize=18, pad=8, fontweight="bold")
    ax.text(
        0.02, 0.96, panel_label,
        transform=ax.transAxes,
        ha='left', va='top',
        fontsize=18, fontweight='bold'
    )

    remove_map_outline(ax)
    return m



def rolling_trend_np(y, window=1000):
    y = np.asarray(y).astype(float)
    N = len(y)
    n_slopes = N - window + 1
    slopes = np.empty(n_slopes)
    years = np.arange(N)[window-1:] 
    for i in range(n_slopes):
        y_win = y[i:i+window]
        t_win = np.arange(window)
        # np.polyfit deg=1 retourne [pente, intercept]
        slope, intercept = np.polyfit(t_win, y_win, 1)
        slopes[i] = slope  # pente par unité de temps (année)
    return slopes, years

def rolling_mean_np(y, window=1500):
    moyennes = []
    for i in range(len(y) - window + 1):
        moyenne = sum(y[i:i+window]) / window
        moyennes.append(moyenne)
    return moyennes


def linear_fit_hac_1d(y, conf=90, nlags=None):
    y = np.asarray(y, dtype=float)
    valid = np.isfinite(y)
    if valid.sum() < 3:
        return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
    yv = y[valid]
    t = np.arange(len(y))[valid]
    X = sm.add_constant(t)
    if nlags is None:
        n = len(yv)
        nlags = int(np.floor(4 * (n / 100)**(2/9)))
    fit = sm.OLS(yv, X).fit(cov_type="HAC", cov_kwds={"maxlags": nlags})
    intercept = fit.params[0]
    slope = fit.params[1]
    stderr_slope = fit.bse[1]
    pvalue = fit.pvalues[1]
    alpha = 1 - conf / 100
    z = norm.ppf(1 - alpha / 2)
    ci_low = slope - z * stderr_slope
    ci_high = slope + z * stderr_slope
    return slope, intercept, stderr_slope, ci_low, ci_high, pvalue


def fit_map_hac(da, conf=90, nlags=None):
    ny = da.sizes["y"]
    nx = da.sizes["x"]
    slope = np.full((ny, nx), np.nan, dtype=float)
    intercept = np.full((ny, nx), np.nan, dtype=float)
    stderr_slope = np.full((ny, nx), np.nan, dtype=float)
    ci_low = np.full((ny, nx), np.nan, dtype=float)
    ci_high = np.full((ny, nx), np.nan, dtype=float)
    pvalue = np.full((ny, nx), np.nan, dtype=float)
    values = da.values  # shape = (time, y, x)
    for j in range(ny):
        for i in range(nx):
            y = values[:, j, i]

            s, a, se, low, high, p = linear_fit_hac_1d(
                y, conf=conf, nlags=nlags
            )

            slope[j, i] = s
            intercept[j, i] = a
            stderr_slope[j, i] = se
            ci_low[j, i] = low
            ci_high[j, i] = high
            pvalue[j, i] = p
    ds_out = xr.Dataset(
        data_vars=dict(
            slope=(("y", "x"), slope),
            intercept=(("y", "x"), intercept),
            stderr_slope=(("y", "x"), stderr_slope),
            ci_low=(("y", "x"), ci_low),
            ci_high=(("y", "x"), ci_high),
            pvalue=(("y", "x"), pvalue),
        ),
        coords={
            "y": da["y"],
            "x": da["x"],
        }
    )
    return ds_out
