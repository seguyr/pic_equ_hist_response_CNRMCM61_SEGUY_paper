"""
Figure 7: Historical spatial response comparison by layer
=========================================================
This script plots, for each depth layer:
- Hist_dd+1000 gain map,
- Hist_dd+3000 - Hist_dd+1000 gain difference map.
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
    load_ohc_2d_layer,
    time_matching,
    anomalies,
    gain,
    boot,
    boot_diff,
    plot_panel,
    make_cmap_norm,
)

# -----------------------------------------------------------------------------
# Parameters
# -----------------------------------------------------------------------------
ECHELLE_OHC = 1e9  # J m^-2 -> GJ m^-2
UNIT = "GJ"

layers = ["0-300", "300-2000", "2000-bottom"]

layer_labels = {
    "0-300": "0–300 m",
    "300-2000": "300–2000 m",
    "2000-bottom": "2000 m–bottom",
}

col_titles = [
    "Hist_dd+1000",
    "Hist_dd+3000 - Hist_dd+1000",
]

row_limits = {
    "0-300": (-2, 2),
    "300-2000": (-3, 3),
    "2000-bottom": (-1.5, 1.5),
}

# -----------------------------------------------------------------------------
# Figure
# -----------------------------------------------------------------------------

proj = ccrs.Robinson(central_longitude=0)

fig = plt.figure(figsize=(16, 15))
gs = fig.add_gridspec(
    nrows=3,
    ncols=3,
    width_ratios=[1, 1, 0.02],
    wspace=0.04,
    hspace=0.18,
)

axes_grid = [[None, None] for _ in range(len(layers))]

for i, layer in enumerate(layers):
    # -------------------------------------------------------------
    # Load one layer
    # -------------------------------------------------------------
    hist_1000 = load_ohc_2d_layer("hist_tot", layer=layer) / ECHELLE_OHC
    hist_3000 = load_ohc_2d_layer("hist_3000", layer=layer) / ECHELLE_OHC
    pic_1000 = load_ohc_2d_layer("pic_tot", layer=layer) / ECHELLE_OHC
    pic_3000 = load_ohc_2d_layer("pic_3000", layer=layer) / ECHELLE_OHC

    # -------------------------------------------------------------
    # Dedrift + anomalies + gain
    # -------------------------------------------------------------
    ohc_dd_1000 = time_matching(hist_1000, pic_1000)
    ohc_dd_3000 = time_matching(hist_3000, pic_3000)

    ohc_dd_1000 = anomalies(ohc_dd_1000)
    ohc_dd_3000 = anomalies(ohc_dd_3000)

    m_ohc_dd_1000 = gain(ohc_dd_1000)
    m_ohc_dd_3000 = gain(ohc_dd_3000)

    # -------------------------------------------------------------
    # Bootstrap on gain maps
    # -------------------------------------------------------------
    hist_gain_1000 = boot(m_ohc_dd_1000)
    diff_hist_gain = boot_diff(m_ohc_dd_1000, m_ohc_dd_3000)

    # -------------------------------------------------------------
    # Subplots
    # -------------------------------------------------------------
    axL = fig.add_subplot(gs[i, 0], projection=proj)
    axR = fig.add_subplot(gs[i, 1], projection=proj)
    cax = fig.add_subplot(gs[i, 2])

    axes_grid[i][0] = axL
    axes_grid[i][1] = axR

    vmin, vmax = row_limits[layer]
    cmap_row, norm_row = make_cmap_norm(vmin, vmax)

    lab_left = f"{chr(97 + 2*i)})"
    lab_right = f"{chr(97 + 2*i + 1)})"

    mL = plot_panel(axL, hist_gain_1000, "", lab_left, cmap_row, norm_row)
    mR = plot_panel(axR, diff_hist_gain, "", lab_right, cmap_row, norm_row)

    cbar = fig.colorbar(
        mR,
        cax=cax,
        orientation="vertical",
        extend="both",
        extendfrac="auto",
    )
    cbar.set_label(f"OHC gain ({UNIT} m$^{-2}$)", fontsize=22, labelpad=16)
    cbar.ax.tick_params(labelsize=20, length=6, width=1.1)


# -----------------------------------------------------------------------------
# Column titles
# -----------------------------------------------------------------------------
p0L = axes_grid[0][0].get_position()
p0R = axes_grid[0][1].get_position()

x_col1 = 0.5 * (p0L.x0 + p0L.x1)
x_col2 = 0.5 * (p0R.x0 + p0R.x1)
y_col_titles = max(p0L.y1, p0R.y1) + 0.02

fig.text(
    x_col1, y_col_titles, col_titles[0],
    ha="center", va="bottom",
    fontsize=20, fontweight="bold"
)
fig.text(
    x_col2, y_col_titles, col_titles[1],
    ha="center", va="bottom",
    fontsize=20, fontweight="bold"
)

# -----------------------------------------------------------------------------
# Row labels
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

plt.savefig(FIG_DIR / "figure6_layers.pdf", bbox_inches="tight")
plt.show()
