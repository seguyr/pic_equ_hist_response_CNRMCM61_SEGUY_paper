"""
Figure 6: Historical spatial response comparison by layer
=========================================================
This script plots, for each depth layer:
- Hist_dd+1000 gain map,
- Hist_dd+3000 - Hist_dd+1000 gain difference map.
Layers:
- 0-300 m
- 300-2000 m
- 2000-bottom
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
    load_area_ocean,
    load_ohc_2d,
    time_matching,
    anomalies,
    gain,
    boot,
    boot_diff,
    plot_panel,
)

# -----------------------------------------------------------------------------
# Parameters
# -----------------------------------------------------------------------------
ECHELLE_OHC = 1e9  # J m^-2 -> GJ m^-2
UNIT = "GJ m$^{-2}$"

LAYERS = ["0-300", "300-2000", "2000-bottom"]

LAYER_LABELS = {
    "0-300": "0–300 m",
    "300-2000": "300–2000 m",
    "2000-bottom": "2000 m–bottom",
}

COL_TITLES = [
    "Hist_dd+1000",
    "Hist_dd+3000 − Hist_dd+1000",
]

ROW_LIMITS = {
    "0-300": (-2, 2),
    "300-2000": (-3, 3),
    "2000-bottom": (-1.5, 1.5),
}

WHITE_FRAC = 0.01


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
def make_cmap_norm(vmin, vmax, white_frac=0.01, n=256):
    """Create a diverging colormap with a narrow white band around zero."""
    white_width = white_frac * (vmax - vmin)

    bounds = np.linspace(vmin, vmax, n)
    colors = plt.cm.RdBu_r(np.linspace(0, 1, n))

    i0_low = np.argmin(np.abs(bounds + white_width))
    i0_high = np.argmin(np.abs(bounds - white_width))
    colors[i0_low:i0_high] = [1, 1, 1, 1]

    cmap = mcolors.ListedColormap(colors)
    norm = mcolors.TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)
    return cmap, norm


def main():
    """Generate the multi-layer spatial gain figure."""

    area_oce = load_area_ocean()
    proj = ccrs.Robinson(central_longitude=0)

    fig = plt.figure(figsize=(16, 15))
    gs = fig.add_gridspec(
        nrows=3,
        ncols=3,
        width_ratios=[1, 1, 0.02],
        wspace=0.04,
        hspace=0.18,
    )

    axes_grid = [[None, None] for _ in range(len(LAYERS))]

    for i, layer in enumerate(LAYERS):
        # -------------------------------------------------------------
        # Load data for this layer
        # -------------------------------------------------------------
        hist_1000 = load_ohc_2d("hist_tot", layer=layer) / ECHELLE_OHC
        hist_3000 = load_ohc_2d("hist_3000", layer=layer) / ECHELLE_OHC
        pic_1000 = load_ohc_2d("pic_tot", layer=layer) / ECHELLE_OHC
        pic_3000 = load_ohc_2d("pic_3000", layer=layer) / ECHELLE_OHC

        # -------------------------------------------------------------
        # Dedrift + anomalies + gain
        # -------------------------------------------------------------
        ohc_dd_1000 = anomalies(time_matching(hist_1000, pic_1000))
        ohc_dd_3000 = anomalies(time_matching(hist_3000, pic_3000))

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

        vmin, vmax = ROW_LIMITS[layer]
        cmap_row, norm_row = make_cmap_norm(vmin, vmax, white_frac=WHITE_FRAC)

        lab_left = f"{chr(97 + 2 * i)})"
        lab_right = f"{chr(97 + 2 * i + 1)})"

        mL = plot_panel(
            axL,
            hist_gain_1000,
            lab_left,
            lab_left,
            area_oce,
            cmap_row,
            norm_row,
        )
        axL.set_title("Hist_dd+1000", fontsize=18, fontweight="bold")

        mR = plot_panel(
            axR,
            diff_hist_gain,
            lab_right,
            lab_right,
            area_oce,
            cmap_row,
            norm_row,
        )
        axR.set_title("Hist_dd+3000 − Hist_dd+1000", fontsize=18, fontweight="bold")

        cbar = fig.colorbar(
            mR,
            cax=cax,
            orientation="vertical",
            extend="both",
            extendfrac="auto",
        )
        cbar.set_label(f"OHC gain ({UNIT})", fontsize=20, labelpad=14)
        cbar.ax.tick_params(labelsize=16, length=6, width=1.0)

    # -------------------------------------------------------------
    # Column titles
    # -------------------------------------------------------------
    p0L = axes_grid[0][0].get_position()
    p0R = axes_grid[0][1].get_position()

    x_col1 = 0.5 * (p0L.x0 + p0L.x1)
    x_col2 = 0.5 * (p0R.x0 + p0R.x1)
    y_col_titles = max(p0L.y1, p0R.y1) + 0.02

    fig.text(
        x_col1,
        y_col_titles,
        COL_TITLES[0],
        ha="center",
        va="bottom",
        fontsize=20,
        fontweight="bold",
    )
    fig.text(
        x_col2,
        y_col_titles,
        COL_TITLES[1],
        ha="center",
        va="bottom",
        fontsize=20,
        fontweight="bold",
    )

    # -------------------------------------------------------------
    # Row labels
    # -------------------------------------------------------------
    for i, layer in enumerate(LAYERS):
        pL = axes_grid[i][0].get_position()
        pR = axes_grid[i][1].get_position()

        y_center = 0.5 * (pL.y0 + pL.y1)
        x_left = min(pL.x0, pR.x0) - 0.035

        fig.text(
            x_left,
            y_center,
            LAYER_LABELS[layer],
            rotation=90,
            ha="center",
            va="center",
            fontsize=18,
            fontweight="bold",
        )

    plt.tight_layout()
    plt.savefig(FIG_DIR / "figure_7.pdf", bbox_inches="tight")
    plt.show()


if __name__ == "__main__":
    main()
