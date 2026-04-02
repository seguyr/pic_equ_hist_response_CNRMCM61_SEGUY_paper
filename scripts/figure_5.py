"""
Figure 5: Historical time series response comparison
====================================================

This script plots:
- pic_tot time series in anomalies relative to 1850-1900,
- pic_3000 time series in anomalies relative to 1850-1900,
- hist_tot time series in anomalies relative to 1850-1900,
- hist_3000 time series in anomalies relative to 1850-1900,
- the difference between historical OHC responses after dedrifting by time matching.
"""

from pathlib import Path
import sys
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from matplotlib.patches import Patch

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
    time_matching,
    anomalies,
    gain,
    get_stats,
    get_scalar_stats,
    COLORS
)

# -----------------------------------------------------------------------------
# Physical constants
# -----------------------------------------------------------------------------
ECHELLE_OHC = 1e21  # convert J to ZJ
UNIT = "ZJ"

REF_START = 1850
REF_END = 1899
P_START = 1995
P_END = 2014
PERIOD_GAIN_START = 1995
PERIOD_GAIN_END = 2014

TIME_HIST = np.arange(165) + REF_START
x_gb, y_gb = 0.7, 0.84


pic_1000 = xr.open_dataset(INTERMEDIATE_DIR / "pic_1000.nc").__xarray_dataarray_variable__
pic_3000 = xr.open_dataset(INTERMEDIATE_DIR / "pic_3000.nc").__xarray_dataarray_variable__
hist_1000 = xr.open_dataset(INTERMEDIATE_DIR / "hist_1000.nc").__xarray_dataarray_variable__
hist_3000 = xr.open_dataset(INTERMEDIATE_DIR / "hist_3000.nc").__xarray_dataarray_variable__

diff_hist_cor = xr.open_dataset(INTERMEDIATE_DIR / "diff_hist_cor.nc").__xarray_dataarray_variable__
m_diff_hist_cor = xr.open_dataset(INTERMEDIATE_DIR / "m_diff_hist_cor.nc").__xarray_dataarray_variable__

rm_diff_hist_cor = m_diff_hist_cor.values[1]
rm_qinf_diff_hist_cor = m_diff_hist_cor.values[0]
rm_qsup_diff_hist_cor = m_diff_hist_cor.values[2]

y = y_gb
x = x_gb + 0.3

# --- Création d'une figure avec deux subplots superposés ---
fig, axes = plt.subplots(2, 1, figsize=(15, 15), sharex=True)
ax1, ax2 = axes
bar_colors = [COLORS["orange_dark"], COLORS["teal_dark"], COLORS["red_dark"], COLORS["blue_dark"], COLORS["red_light"]]
n_bars = 5
bar_labels = ['Hist_dd+1000', 'Hist+3000', 'Pic+1000', 'PiC+3000', "Hist_dd+3000 - Hist_dd+1000"]
panel_labels = ["a)", "b)"]

# ============================================================
# --- SUBPLOT 1 : piControl & historique bruts ---
# ============================================================
time = pic_1000.time + REF_START
ax1.plot(time, pic_1000.sel(stats='mean').values, color=COLORS["red_dark"])
ax1.fill_between(time, pic_1000.sel(stats='lower').values, pic_1000.sel(stats='upper').values, 
                 color=COLORS["red_light"], alpha=0.8)

ax1.plot(time, pic_3000.sel(stats='mean').values, color=COLORS["blue_dark"])
ax1.fill_between(time, pic_3000.sel(stats='lower').values, pic_3000.sel(stats='upper').values, 
                 color=COLORS["blue_light"], alpha=0.8)

ax1.plot(time, hist_1000.sel(stats='mean').values, color=COLORS["orange_dark"])
ax1.fill_between(time, hist_1000.sel(stats='lower').values, hist_1000.sel(stats='upper').values , 
                 color=COLORS["orange_light"], alpha=0.8)

ax1.plot(time, hist_3000.sel(stats='mean').values, color=COLORS["teal_dark"])
ax1.fill_between(time, hist_3000.sel(stats='lower').values, hist_3000.sel(stats='upper').values, 
                 color=COLORS["teal_light"], alpha=0.8)

ax1.set_ylabel(f"OHC / {REF_START}-{REF_END} ({UNIT})", fontsize=20)
ax1.set_title("CNRM-CM6.1 working ensembles", fontsize=25, fontweight="bold")
ax1.tick_params(axis='both', labelsize=20)
ax1.axvspan(REF_START, REF_END, color='grey', alpha=0.2)

ax1.grid()

# --- label a) (hors cadre, à gauche)
ax1.text(
    -0.08, 1.02, panel_labels[0],
    transform=ax1.transAxes,
    ha="right", va="bottom",
    fontsize=22,
    fontweight="bold",
    clip_on=False
)

# ============================================================
# --- SUBPLOT 2 : Différence hist_corr_anom ---
# ============================================================
print("Différence hist_corr_anom soustraction ==> global")
time = diff_hist_cor.time + REF_START
b_diff = diff_hist_cor
r_diff = rm_diff_hist_cor
qinf_diff = rm_qinf_diff_hist_cor
qsup_diff = rm_qsup_diff_hist_cor

ax2.plot(time, b_diff.sel(stats='mean').values, color="grey")
ax2.fill_between(time, b_diff.sel(stats='lower').values, b_diff.sel(stats='upper').values, 
                 color='red', alpha=0.3)
ax2.hlines(r_diff, 1995, 2014, colors='purple', linewidth=3)
ax2.errorbar(x=2005, y=r_diff, 
             yerr=[[r_diff - qinf_diff], [qsup_diff - r_diff]], 
             fmt='o', color='purple', capsize=5)
ax2.set_ylabel(f"OHC / {REF_START}-{REF_END} ({UNIT})", fontsize=20)
ax2.set_xlabel("Year ", fontsize=20)
ax2.tick_params(axis='both', labelsize=20)
ax2.axvspan(REF_START, REF_END, color='grey', alpha=0.2)
ax2.axvspan(P_START, P_END, color='grey', alpha=0.2)
ax2.grid()
ax2.axhline(0, color='green', linewidth=3)
ax2.set_title("Hist_dd+3000 - Hist_dd+1000", fontsize=25, fontweight="bold")

# --- label b) (hors cadre, à gauche)
ax2.text(
    -0.08, 1.02, panel_labels[1],
    transform=ax2.transAxes,
    ha="right", va="bottom",
    fontsize=22,
    fontweight="bold",
    clip_on=False
)
bar_handles = [
    Patch(facecolor=bar_colors[k], edgecolor="none", label=bar_labels[k])
    for k in range(n_bars)
]

fig.legend(
    handles=bar_handles,
    loc="upper center",
    ncol=5,
    frameon=False,
    fontsize=20,
    bbox_to_anchor=(0.5, -0.01)
)

plt.tight_layout(h_pad=3)
plt.savefig(FIG_DIR / "figure_5.pdf", bbox_inches="tight")

