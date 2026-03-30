
"""
Figure 6: Historical spatial response comparison
====================================================

This script plots:
- hist_tot map gain,
- hist_3000 map gain,
- the gain difference between historical OHC spatial responses after dedrifting by time matching.
"""


import numpy as np
import xarray as xr
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from dask.diagnostics import ProgressBar
import datetime 
from scipy import stats
from matplotlib.cm import get_cmap
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors


# -----------------------------------------------------------------------------
# Make the project root importable
# -----------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.append(str(PROJECT_ROOT))

# -----------------------------------------------------------------------------
# Output directory for figures
# -----------------------------------------------------------------------------
FIG_DIR = PROJECT_ROOT / "figures"
FIG_DIR.mkdir(exist_ok=True)

# -----------------------------------------------------------------------------
# Imports from project utilities
# -----------------------------------------------------------------------------
from functions.utils import (
    load_ohc_2d_ensembles,
    time_matching,
    anomalies,
    gain,
    boot,
    boot_diff,
    get_stats,
    get_scalar_stats,
    COLORS
)

# -----------------------------------------------------------------------------
# Physical constants
# -----------------------------------------------------------------------------
ECHELLE_OHC = 1e9  # convert J to GJ
UNIT = "GJ"

REF_START = 1850
REF_END = 1899

PERIOD_GAIN_START = 1995
PERIOD_GAIN_END = 2014

TIME_HIST = np.arange(165) + REF_START

# -----------------------------------------------------------------------------
# Calcul
# -----------------------------------------------------------------------------

hist_1000 = load_ohc_2d_ensembles("hist_tot") / ECHELLE_OHC
hist_3000 = load_ohc_2d_ensembles("hist_3000") / ECHELLE_OHC
pic_1000 = load_ohc_2d_ensembles("pic_tot") / ECHELLE_OHC
pic_3000 = load_ohc_2d_ensembles("pic_3000") / ECHELLE_OHC

# Dedrift by time matching
OHC_dd_1000 = time_matching(hist_1000, pic_1000)
OHC_dd_3000 = time_matching(hist_3000, pic_3000)

# Anomalies relative to first 50 years
OHC_dd_1000 = anomalies(OHC_dd_1000)
OHC_dd_3000 = anomalies(OHC_dd_3000)

# Mean gain over final 20 years
m_OHC_dd_1000 = gain(OHC_dd_1000)
m_OHC_dd_3000 = gain(OHC_dd_3000)

# Bootstrap
hist_cor_anom_1000 = boot(OHC_dd_1000)
hist_cor_anom_3000 = boot(OHC_dd_3000)
diff_hist_cor = boot_diff(OHC_dd_1000, OHC_dd_3000)

pic_1000_low, pic_1000_mean, pic_1000_up = get_stats(pic_1000)
pic_3000_low, pic_3000_mean, pic_3000_up = get_stats(pic_3000)

hist_1000_low, hist_1000_mean, hist_1000_up = get_stats(hist_1000)
hist_3000_low, hist_3000_mean, hist_3000_up = get_stats(hist_3000)

diff_low, diff_mean, diff_up = get_stats(diff_hist_cor)



vmin, vmax = -5, 5
white_width = 0.01 * (vmax - vmin)

bounds = np.linspace(vmin, vmax, 256)
colors = plt.cm.RdBu_r(np.linspace(0, 1, 256))
i0_low  = np.argmin(np.abs(bounds + white_width))
i0_high = np.argmin(np.abs(bounds - white_width))
colors[i0_low:i0_high] = [1, 1, 1, 1]

cmap_custom = mcolors.ListedColormap(colors)
norm = mcolors.TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)


def get_stats(ds):
    da_low  = ds.sel(stats="lower").__xarray_dataarray_variable__
    da_mean = ds.sel(stats="mean").__xarray_dataarray_variable__
    da_up   = ds.sel(stats="upper").__xarray_dataarray_variable__
    return da_low, da_mean, da_up

def remove_map_outline(ax):
    """Enlève complètement le contour noir (outline) de la projection."""
    if hasattr(ax, "outline_patch"):
        ax.outline_patch.set_visible(False)
    if "geo" in ax.spines:
        ax.spines["geo"].set_visible(False)
    ax.patch.set_edgecolor("none")
    ax.patch.set_linewidth(0)

def plot_panel(ax, ds, title, label, lon2d, lat2d):
    low, mean, up = get_stats(ds)

    ax.set_global()

    # Coastlines/borders propres (publi)
    ax.coastlines(linewidth=0.55)
    ax.add_feature(cfeature.BORDERS, linewidth=0.35, linestyle=':')

    # Champ principal
    cf = ax.pcolormesh(
        lon2d, lat2d, mean,
        cmap=cmap_custom, norm=norm,
        transform=ccrs.PlateCarree()
    )

    # Hachures : 0 ∈ IC (non significatif)
    mask_non_sig = np.ma.masked_where((low >= 0) | (up <= 0), mean)
    ax.contourf(
        lon2d, lat2d, mask_non_sig,
        hatches=['///'],
        colors='none',
        transform=ccrs.PlateCarree(),
        zorder=2
    )
    
    hatch = ax.contourf(
    lon2d, lat2d, mask_non_sig,
    hatches=['///'],
    colors='none',
    transform=ccrs.PlateCarree(),
    zorder=2
    )

    # rendre les hachures grises
    for coll in hatch.collections:
        coll.set_edgecolor('lightgray')
        coll.set_linewidth(0.0)

    # Titres + labels a) b) c)
    ax.set_title(title, fontsize=16, pad=8, fontweight="bold")
    ax.text(
        0.02, 0.96, label,
        transform=ax.transAxes,
        ha='left', va='top',
        fontsize=18, fontweight='bold'
    )

    remove_map_outline(ax)

    return cf

# -----------------------------
# FIGURE 
# -----------------------------
proj = ccrs.Robinson(central_longitude=0)


fig = plt.figure(figsize=(12.5, 16))
gs = fig.add_gridspec(
    nrows=3, ncols=2,
    width_ratios=[1, 0.03],
    wspace=0.03, hspace=0.12
)

axes = [fig.add_subplot(gs[i, 0], projection=proj) for i in range(3)]
cax  = fig.add_subplot(gs[:, 1])

cf = plot_panel(axes[0], boot_1000, "Hist_dd +1000", "a)", lon2d_oce, lat2d_oce)
plot_panel(axes[1], boot_3000, "Hist_dd +3000", "b)", lon2d_oce, lat2d_oce)
plot_panel(axes[2], boot_diff, "Hist_dd +3000 − Hist_dd +1000", "c)", lon2d_oce, lat2d_oce)


cbar = fig.colorbar(
    cf, cax=cax,
    orientation='vertical',
    extend='both',
    extendfrac='auto'
)

cbar.set_label("OHC gain (GJ m$^{-2}$)", fontsize=24, labelpad=20)
cbar.ax.tick_params(labelsize=22, length=7, width=1.2)

plt.show()
