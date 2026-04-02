"""
Figure 1: Experiment protocol
=============================

This script plots:
- the reference piControl global ocean heat content (OHC),
- the reference AMOC time series,
- the historical ensemble members branched after 1000 years spin-up,
- the historical ensemble members branched after 3000 years spin-up.

"""

from pathlib import Path
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
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
    load_branching_years,
)


# -----------------------------------------------------------------------------
# Physical constants
# -----------------------------------------------------------------------------
ECHELLE_OHC = 1e21  # ZJ

# -----------------------------------------------------------------------------
# Time axes
# -----------------------------------------------------------------------------
TIME_HIST = np.arange(165)   # 0 to 164

# -----------------------------------------------------------------------------
# Shaded branching windows
# -----------------------------------------------------------------------------
YINI_1000 = 1883 - 1850
YEND_1000 = 2831 - 1850
YINI_3000 = 4462 - 1850
YEND_3000 = 4685 - 1850


amoc_pic_gb = xr.open_dataset(INTERMEDIATE_DIR / "amoc_pic_gb.nc").msftyz
ohc_pic_gb = xr.open_dataset(INTERMEDIATE_DIR / "ohc_pic_gb.nc").thetao

ohc_hist_1000 = xr.open_dataset(INTERMEDIATE_DIR / "ohc_hist_1000.nc").__xarray_dataarray_variable__
ohc_hist_3000 = xr.open_dataset(INTERMEDIATE_DIR / "ohc_hist_3000.nc").__xarray_dataarray_variable__
ohc_pic_1000 = xr.open_dataset(INTERMEDIATE_DIR / "ohc_pic_1000.nc").__xarray_dataarray_variable__
ohc_pic_3000 = xr.open_dataset(INTERMEDIATE_DIR / "ohc_pic_3000.nc").__xarray_dataarray_variable__


starty_tot, starty_3000 = load_branching_years()
time = np.linspace(0,2999,3000)
time_hist = np.linspace(0,164,165)

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------
data = ohc_pic_gb
data_tas = amoc_pic_gb
cmap2 = get_cmap('Oranges')
cmap1 = get_cmap('Blues')
#tracer
fig, ax1 = plt.subplots(figsize=(15, 10))
ax1.plot(data.time, data, color="black", label="OHC", linewidth = 3)
cpt = 0

for k in starty_tot:
    ax1.plot(k - 1850, 0, marker = 'o', color = cmap2(cpt/28), markersize = 10, transform = ax1.get_xaxis_transform())
    ohc_hist_1000[cpt]['time'] = time_hist
    hist = ohc_hist_1000.isel(ensemble=cpt)
    ax1.plot(time_hist + k - 1850, hist,'--', color = cmap2(cpt/28),linewidth = 2)
    cpt +=1
cpt = 0
for k in  starty_3000:
    ax1.plot(k - 1850, 0, marker = 'o', color = cmap1(cpt/9), markersize = 10, transform = ax1.get_xaxis_transform())
    ohc_hist_3000[cpt]['time'] = time_hist
    hist2 = ohc_hist_3000.isel(ensemble=cpt)
    ax1.plot(time_hist + k - 1850, hist2,'--', color = cmap1(cpt/9),linewidth = 2)
    cpt +=1

ax1.axvspan(YINI_1000,YEND_1000, color = 'orange', alpha = 0.3, label = "after 1000 years Spin-up")
ax1.axvspan(YINI_3000,YEND_3000, color = 'blue', alpha = 0.3, label = "after 3000 years Spin-up")
ax1.set_ylabel("Global OHC (1e3 ZJ)", color="black", fontsize=30)
ax1.set_xlabel("Time (simulated years)", fontsize=30)
ax1.tick_params(axis='y', labelcolor="black", labelsize=20)
ax1.tick_params(axis='x', labelcolor="black", labelsize=20)
ax1.grid(True)
ax1.set_xlim(0,3000)

ax2 = ax1.twinx()

ax2.plot(
    data_tas.year,
    data_tas,
    color="grey",
    linewidth=2,
    alpha=0.4,
    label="AMOC"
)

ax2.set_ylabel("AMOC (Sv)", fontsize=30, color="grey")
ax2.tick_params(axis='y', labelcolor="grey", labelsize=20)

ax2.set_zorder(0)    
ax1.set_zorder(1)
ax1.patch.set_visible(False)


# --- Légende personnalisée ---
legend_elements = [
    Patch(facecolor='orange', edgecolor='orange', label='hist+1000'), # rectangle rouge
    Patch(facecolor='blue', edgecolor='blue', label='hist+3000'), # rectangle violet
    Line2D([0], [0], marker='o', color='black', label='branching time',
           markersize=10, linestyle='None'), # rond noir
    Line2D([0], [0], color='black', lw=3, label='Reference simulation OHC') # trait noir
]
legend_elements.append(
    Line2D([0], [0], color='grey', lw=2, alpha=0.4, label='Reference simulation AMOC')  # trait gris
)

ax1.legend(handles=legend_elements, loc='center', bbox_to_anchor=(0.5,-0.2),ncol = 2, fontsize=20)
ax1.set_title('Experience Protocol', fontsize = 35)

plt.savefig(FIG_DIR / "figure_1.pdf", bbox_inches="tight")
plt.savefig(FIG_DIR / "figure_1.png", bbox_inches="tight")

