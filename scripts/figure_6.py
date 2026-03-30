
"""
Figure 6: Historical spatial response comparison
====================================================
This script plots:
- hist_tot map gain,
- hist_3000 map gain,
- the gain difference between historical OHC spatial responses after
dedrifting by time matching.
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
TIME_HIST = np.arange(165) + REF_START

# -----------------------------------------------------------------------------
# Visualisation
# -----------------------------------------------------------------------------

VMIN, VMAX = -5, 5
WHITE_WIDTH = 0.01 * (VMAX - VMIN)

bounds = np.linspace(VMIN, VMAX, 256)
colors = plt.cm.RdBu_r(np.linspace(0, 1, 256))
i0_low = np.argmin(np.abs(bounds + WHITE_WIDTH))
i0_high = np.argmin(np.abs(bounds - WHITE_WIDTH))
colors[i0_low:i0_high] = [1, 1, 1, 1]

CMAP_CUSTOM = mcolors.ListedColormap(colors)
NORM = mcolors.TwoSlopeNorm(vmin=VMIN, vcenter=0, vmax=VMAX)

# -----------------------------------------------------------------------------
# Calcul
# -----------------------------------------------------------------------------

hist_1000 = load_ohc_2d("hist_tot", layer="all") / ECHELLE_OHC
hist_3000 = load_ohc_2d("hist_3000", layer="all") / ECHELLE_OHC
pic_1000 = load_ohc_2d("pic_tot", layer="all") / ECHELLE_OHC
pic_3000 = load_ohc_2d("pic_3000", layer="all") / ECHELLE_OHC

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
hist_gain_1000 = boot(m_ohc_dd_1000)
hist_gain_3000 = boot(m_ohc_dd_3000)
diff_hist_gain = boot_diff(m_ohc_dd_1000, m_ohc_dd_3000)

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

cf = plot_panel(axes[0], hist_gain_1000, "Hist_dd +1000", "a)", CMAP_CUSTOM, NORM)
plot_panel(axes[1], hist_gain_3000, "Hist_dd +3000", "b)", CMAP_CUSTOM, NORM)
plot_panel(axes[2], diff_hist_gain, "Hist_dd +3000 − Hist_dd +1000", "c)", CMAP_CUSTOM, NORM)

cbar = fig.colorbar(
    cf, cax=cax,
    orientation='vertical',
    extend='both',
    extendfrac='auto'
)

cbar.set_label(f"OHC gain ({UNIT} m$^{-2}$)", fontsize=24, labelpad=20)
cbar.ax.tick_params(labelsize=22, length=7, width=1.2)
plt.savefig(FIG_DIR / "figure_6.pdf", bbox_inches="tight")

