
"""
Figure 6: Historical spatial response comparison
====================================================

This script plots:
- hist_tot map gain,
- hist_3000 map gain,
- the gain difference between historical OHC spatial responses after dedrifting by time matching.
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
    plot_panel,
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

vmin, vmax = -5, 5
white_width = 0.01 * (vmax - vmin)

bounds = np.linspace(vmin, vmax, 256)
colors = plt.cm.RdBu_r(np.linspace(0, 1, 256))
i0_low  = np.argmin(np.abs(bounds + white_width))
i0_high = np.argmin(np.abs(bounds - white_width))
colors[i0_low:i0_high] = [1, 1, 1, 1]

cmap_custom = mcolors.ListedColormap(colors)
norm = mcolors.TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)



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

cf = plot_panel(axes[0], hist_cor_anom_1000, "Hist_dd +1000", "a)")
plot_panel(axes[1], hist_cor_anom_3000, "Hist_dd +3000", "b)")
plot_panel(axes[2], diff_hist_cor, "Hist_dd +3000 − Hist_dd +1000", "c)")

cbar = fig.colorbar(
    cf, cax=cax,
    orientation='vertical',
    extend='both',
    extendfrac='auto'
)

cbar.set_label("OHC gain (GJ m$^{-2}$)", fontsize=24, labelpad=20)
cbar.ax.tick_params(labelsize=22, length=7, width=1.2)

plt.show()
