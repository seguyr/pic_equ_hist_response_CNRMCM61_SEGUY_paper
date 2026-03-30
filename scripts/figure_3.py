"""
Figure 3: piControl OHC drift maps
==================================
This script plots:
- the OHC trend map over the first 1000 simulated years,
- the OHC trend map over the last 400 simulated years.
Trends are estimated gridpoint by gridpoint using a linear fit with
HAC-corrected confidence intervals.
"""

from pathlib import Path
import sys
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cartopy.crs as ccrs

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
from functions.utils import (  # noqa: E402
    DIR_PIC,
    load_area_ocean,
    fit_map_hac,
    plot_panel,
    STATS_COORD,
)

# -----------------------------------------------------------------------------
# Plot style
# -----------------------------------------------------------------------------
plt.rc("font", size=12.5)
plt.rcParams["figure.dpi"] = 300

# -----------------------------------------------------------------------------
# Parameters
# -----------------------------------------------------------------------------
ECHELLE = 1e9  # J m^-2 -> GJ m^-2
UNIT = "GJ m$^{-2}$ century$^{-1}$"
CONF = 90

OHC_PIC_FILE = DIR_PIC / "pic_3000.nc"

# -----------------------------------------------------------------------------
# Load data
# -----------------------------------------------------------------------------
ds_pic = xr.open_dataset(OHC_PIC_FILE)
varname = list(ds_pic.data_vars)[0]
ohc_pic = ds_pic[varname]

# harmonize time coordinate if needed
if "time" in ohc_pic.dims and ohc_pic.sizes["time"] == 3000:
    ohc_pic = ohc_pic.assign_coords(time=np.arange(3000))

# convert to GJ m^-2
ohc_pic = ohc_pic / ECHELLE

# first 1000 years and last 400 years
ohc_1000 = ohc_pic.isel(time=slice(0, 1000))
ohc_400 = ohc_pic.isel(time=slice(2600, 3000))

# -----------------------------------------------------------------------------
# HAC fits
# -----------------------------------------------------------------------------
fit_map_1000 = fit_map_hac(ohc_1000, conf=CONF)
fit_map_400 = fit_map_hac(ohc_400, conf=CONF)

# convert slopes from per year to per century
slope_1000 = fit_map_1000["slope"] * 100
ci_low_1000 = fit_map_1000["ci_low"] * 100
ci_high_1000 = fit_map_1000["ci_high"] * 100

slope_400 = fit_map_400["slope"] * 100
ci_low_400 = fit_map_400["ci_low"] * 100
ci_high_400 = fit_map_400["ci_high"] * 100

# -----------------------------------------------------------------------------
# Convert HAC output to the format expected by plot_panel()
# -----------------------------------------------------------------------------
stats_1000 = xr.concat(
    [
        ci_low_1000,
        slope_1000,
        ci_high_1000,
        xr.full_like(slope_1000, np.nan),
        xr.full_like(slope_1000, np.nan),
    ],
    dim=xr.IndexVariable("stats", STATS_COORD),
)

stats_400 = xr.concat(
    [
        ci_low_400,
        slope_400,
        ci_high_400,
        xr.full_like(slope_400, np.nan),
        xr.full_like(slope_400, np.nan),
    ],
    dim=xr.IndexVariable("stats", STATS_COORD),
)

# -----------------------------------------------------------------------------
# Colormap
# -----------------------------------------------------------------------------

VMIN, VMAX = -3, 3
WHITE_WIDTH = 0.01 * (VMAX - VMIN)

bounds = np.linspace(VMIN, VMAX, 256)
colors = plt.cm.RdBu_r(np.linspace(0, 1, 256))
i0_low = np.argmin(np.abs(bounds + WHITE_WIDTH))
i0_high = np.argmin(np.abs(bounds - WHITE_WIDTH))
colors[i0_low:i0_high] = [1, 1, 1, 1]

CMAP_CUSTOM = mcolors.ListedColormap(colors)
NORM = mcolors.TwoSlopeNorm(vmin=VMIN, vcenter=0, vmax=VMAX)

# -----------------------------------------------------------------------------
# Figure
# -----------------------------------------------------------------------------
proj = ccrs.Robinson(central_longitude=0)

fig = plt.figure(figsize=(15, 5))
gs = fig.add_gridspec(
    nrows=2,
    ncols=2,
    height_ratios=[1, 0.05],
    hspace=0.1,
    wspace=0.03,
)

axes = [
    fig.add_subplot(gs[0, 0], projection=proj),
    fig.add_subplot(gs[0, 1], projection=proj),
]

cax = fig.add_subplot(gs[1, :])

cf = plot_panel(axes[0], stats_1000, "Pic+1000", "a)", CMAP_CUSTOM, NORM)
plot_panel(axes[1], stats_400, "Pic+3000", "b)", CMAP_CUSTOM, NORM)

cbar = fig.colorbar(
    cf,
    cax=cax,
    orientation="horizontal",
    extend="both",
    extendfrac="auto",
)

cbar.ax.xaxis.set_label_position("bottom")
cbar.ax.xaxis.set_ticks_position("bottom")
cbar.set_label(f"OHC trend ({UNIT})", fontsize=20, labelpad=20)
cbar.ax.tick_params(labelsize=20, length=7, width=1.2)

plt.tight_layout()
plt.savefig(FIG_DIR / "figure_3.pdf", bbox_inches="tight")

plt.show()
