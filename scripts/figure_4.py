"""
Figure 4 : piControl OHC drift maps by layer
========================================
This script plots, for each depth layer:
- OHC trend over the first 1000 simulated years,
- OHC trend over the last 400 simulated years.
Trends are estimated gridpoint by gridpoint using a linear fit with
HAC-corrected confidence intervals.
"""

from pathlib import Path
import sys
import numpy as np
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
from functions.utils import (
    load_ohc_pic_global,
    fit_map_hac,
    plot_panel_hac,
)

# -----------------------------------------------------------------------------
# Parameters
# -----------------------------------------------------------------------------
ECHELLE = 1e9  # J m^-2 -> GJ m^-2
CONF = 90

layers = ["0_300", "300_2000", "2000_btm"]

layer_labels = {
    "0_300": "0–300 m",
    "300_2000": "300–2000 m",
    "2000_btm": "2000 m–bottom",
}

col_titles = [
    "PiC+1000",
    "PiC+3000",
]

row_limits = {
    "0_300": (-0.5, 0.5),
    "300_2000": (-1, 1),
    "2000_btm": (-1.5, 1.5),
}


# -----------------------------------------------------------------------------
# Calcul des stats par couche
# -----------------------------------------------------------------------------
results = {}

for layer in layers:
    ohc_pic = load_ohc_pic_global(layer) / ECHELLE
    ohc_pic = ohc_pic.assign_coords(time=np.arange(ohc_pic.sizes["time"]))

    ohc_1000 = ohc_pic.isel(time=slice(0, 1000))
    ohc_400 = ohc_pic.isel(time=slice(2600, 3000))

    fit_map_1000 = fit_map_hac(ohc_1000, conf=CONF)
    fit_map_400 = fit_map_hac(ohc_400, conf=CONF)

    results[layer] = {
        "slope_1000": fit_map_1000["slope"] * 100,
        "ci_low_1000": fit_map_1000["ci_low"] * 100,
        "ci_high_1000": fit_map_1000["ci_high"] * 100,
        "slope_400": fit_map_400["slope"] * 100,
        "ci_low_400": fit_map_400["ci_low"] * 100,
        "ci_high_400": fit_map_400["ci_high"] * 100,
    }


# -----------------------------------------------------------------------------
# Figure
# -----------------------------------------------------------------------------
proj = ccrs.Robinson(central_longitude=0)

fig = plt.figure(figsize=(16, 15))
gs = fig.add_gridspec(
    nrows=3, ncols=3,
    width_ratios=[1, 1, 0.025],
    wspace=0.04, hspace=0.18
)

axes_grid = [[None, None] for _ in range(len(layers))]

for i, layer in enumerate(layers):
    axL = fig.add_subplot(gs[i, 0], projection=proj)
    axR = fig.add_subplot(gs[i, 1], projection=proj)
    cax = fig.add_subplot(gs[i, 2])

    axes_grid[i][0] = axL
    axes_grid[i][1] = axR

    vmin, vmax = row_limits[layer]
    cmap_row, norm_row = make_cmap_norm(vmin, vmax)

    lab_left = f"{chr(97 + 2*i)})"
    lab_right = f"{chr(97 + 2*i + 1)})"

    mL = plot_panel_hac(
        axL, col_titles[0], lab_left,
        results[layer]["ci_low_1000"],
        results[layer]["slope_1000"],
        results[layer]["ci_high_1000"],
        cmap_row, norm_row
    )

    mR = plot_panel_hac(
        axR, col_titles[1], lab_right,
        results[layer]["ci_low_400"],
        results[layer]["slope_400"],
        results[layer]["ci_high_400"],
        cmap_row, norm_row
    )

    cbar = fig.colorbar(
        mR, cax=cax,
        orientation="vertical",
        extend="both",
        extendfrac="auto"
    )
    cbar.set_label("OHC trend (GJ m$^{-2}$ century$^{-1}$)", fontsize=15, labelpad=14)
    cbar.ax.tick_params(labelsize=16, length=6, width=1.1)


# -----------------------------------------------------------------------------
# Titres de colonnes
# -----------------------------------------------------------------------------
p0L = axes_grid[0][0].get_position()
p0R = axes_grid[0][1].get_position()

x_col1 = 0.5 * (p0L.x0 + p0L.x1)
x_col2 = 0.5 * (p0R.x0 + p0R.x1)
y_col_titles = max(p0L.y1, p0R.y1) + 0.02

fig.text(x_col1, y_col_titles, col_titles[0],
         ha="center", va="bottom", fontsize=20, fontweight="bold")
fig.text(x_col2, y_col_titles, col_titles[1],
         ha="center", va="bottom", fontsize=20, fontweight="bold")


# -----------------------------------------------------------------------------
# Labels verticaux des couches
# -----------------------------------------------------------------------------
for i, layer in enumerate(layers):
    pL = axes_grid[i][0].get_position()
    pR = axes_grid[i][1].get_position()

    y_center = 0.5 * (pL.y0 + pL.y1)
    x_left = min(pL.x0, pR.x0) - 0.035

    fig.text(
        x_left, y_center,
        layer_labels[layer],
        rotation=90,
        ha="center", va="center",
        fontsize=20, fontweight="bold"
    )

plt.tight_layout()

plt.savefig(FIG_DIR / "figure_4.pdf", bbox_inches="tight")
plt.show()
