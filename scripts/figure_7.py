
"""
Figure 7: Historical spatial response comparison by layers
==========================================================
This script plots, for each layer:
- Hist_dd +1000
- Hist_dd +3000 - Hist_dd +1000
"""

from pathlib import Path
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import xarray as xr

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
INTERMEDIATE_DIR = PROJECT_ROOT / "data/data_plot"

# -----------------------------------------------------------------------------
# Imports from project utilities
# -----------------------------------------------------------------------------
from functions.utils import (
    load_area_ocean,
    nemo_lon_lat,
    remove_map_outline,
    get_stats,
    make_cmap_norm,
)

# -----------------------------------------------------------------------------
# Physical constants
# -----------------------------------------------------------------------------
UNIT = "GJ"

# -----------------------------------------------------------------------------
# Layers and labels
# -----------------------------------------------------------------------------
layers = ["0_300", "300_2000", "2000_btm"]

layer_labels = {
    "0_300": "0–300 m",
    "300_2000": "300–2000 m",
    "2000_btm": "2000 m–bottom",
}

col_titles = [
    "Hist_dd +1000",
    "Hist_dd +3000 − Hist_dd +1000",
]

# Limites par couche
row_limits = {
    "0_300": (-2.0, 2.0),
    "300_2000": (-3.0, 3.0),
    "2000_btm": (-1.5, 1.5),
}

WHITE_FRAC = 0.01


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


# -----------------------------------------------------------------------------
# Load data
# -----------------------------------------------------------------------------
data_1000 = {}
data_diff = {}

for layer in layers:
    data_1000[layer] = xr.open_dataset(
        INTERMEDIATE_DIR / f"boot_1000_{layer}.nc"
    ).__xarray_dataarray_variable__

    data_diff[layer] = xr.open_dataset(
        INTERMEDIATE_DIR / f"boot_diff_{layer}.nc"
    ).__xarray_dataarray_variable__

# -----------------------------------------------------------------------------
# Figure
# -----------------------------------------------------------------------------
proj = ccrs.Robinson(central_longitude=0)

fig = plt.figure(figsize=(15, 15))
gs = fig.add_gridspec(
    nrows=3,
    ncols=3,
    width_ratios=[1, 1, 0.035],
    wspace=0.04,
    hspace=0.14,
)

axes_grid = [[None, None] for _ in range(3)]

for i, layer in enumerate(layers):
    ax_left = fig.add_subplot(gs[i, 0], projection=proj)
    ax_right = fig.add_subplot(gs[i, 1], projection=proj)
    cax = fig.add_subplot(gs[i, 2])

    axes_grid[i][0] = ax_left
    axes_grid[i][1] = ax_right

    vmin, vmax = row_limits[layer]
    cmap_row, norm_row = make_cmap_norm(vmin, vmax, white_frac=WHITE_FRAC)

    label_left = f"{chr(97 + 2 * i)})"
    label_right = f"{chr(97 + 2 * i + 1)})"

    cf_left = plot_panel(
        ax_left,
        data_1000[layer],
        col_titles[0],
        label_left,
        cmap_row,
        norm_row,
    )

    cf_right = plot_panel(
        ax_right,
        data_diff[layer],
        col_titles[1],
        label_right,
        cmap_row,
        norm_row,
    )

    cbar = fig.colorbar(
        cf_right,
        cax=cax,
        orientation="vertical",
        extend="both",
        extendfrac="auto",
    )
    cbar.set_label(f"OHC gain ({UNIT} m$^{{-2}}$)", fontsize=18, labelpad=14)
    cbar.ax.tick_params(labelsize=16, length=6, width=1.0)

# -----------------------------------------------------------------------------
# Row labels on the left
# -----------------------------------------------------------------------------
for i, layer in enumerate(layers):
    pos = axes_grid[i][0].get_position()
    y_center = 0.5 * (pos.y0 + pos.y1)
    x_text = pos.x0 - 0.045

    fig.text(
        x_text,
        y_center,
        layer_labels[layer],
        rotation=90,
        ha="center",
        va="center",
        fontsize=18,
        fontweight="bold",
    )

# -----------------------------------------------------------------------------
# Save
# -----------------------------------------------------------------------------
plt.savefig(FIG_DIR / "figure_7.pdf", bbox_inches="tight")
plt.savefig(FIG_DIR / "figure_7.png", bbox_inches="tight")
plt.show()