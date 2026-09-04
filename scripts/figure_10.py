
"""
Figure 10: Dedrift method impact on results
====================================================

"""

from pathlib import Path
import sys
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch
from scipy import stats
from PIL import Image
import cartopy.crs as ccrs
import cartopy.feature as cfeature

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
INTERMEDIATE_DIR = PROJECT_ROOT / "data/data_plot/fig10"


diff_30_TM = xr.open_mfdataset(INTERMEDIATE_DIR / "diff_30_TM.nc").__xarray_dataarray_variable__
m_diff_30_TM = xr.open_dataset(INTERMEDIATE_DIR / "m_diff_30_TM.nc").__xarray_dataarray_variable__
diff_30_DL = xr.open_mfdataset(INTERMEDIATE_DIR / "diff_30_DL.nc").__xarray_dataarray_variable__
m_diff_30_DL = xr.open_dataset(INTERMEDIATE_DIR / "m_diff_30_DL.nc").__xarray_dataarray_variable__
var = "ohc"
IC = "90%
unit = 'ZJ'
uref_start = 1850
ref_end = 1900

# === Tracé ===
y = 0.7
x = 0.84
time = diff_30_TM.time + 1850
fig, ax = plt.subplots(1, 1, figsize=(15, 10), sharey= False)

ax.plot(time, diff_30_TM.sel(stats='mean').values, color="purple", linewidth=3, label='Time Matching dedrift with 29 members')
ax.fill_between(time, diff_30_TM.sel(stats='lower').values, diff_30_TM.sel(stats='upper').values, color="purple", alpha=0.6)

ax.plot(time, diff_30_DL.sel(stats='mean').values, color="lightgreen", linewidth=3, label='Linear Dedrift dedrift with 29 members')
ax.fill_between(time, diff_30_DL.sel(stats='lower').values, diff_30_DL.sel(stats='upper').values, color="lightgreen", alpha=0.5)

ax.hlines(m_diff_30_TM.sel(stats='mean').values, 1995, 2014, colors='purple', linewidth=3)
ax.errorbar(x=2005, y=m_diff_30_TM.sel(stats='mean').values, 
             yerr=[[m_diff_30_TM.sel(stats='mean').values  - m_diff_30_TM.sel(stats='lower').values ], [m_diff_30_TM.sel(stats='upper').values  - m_diff_30_TM.sel(stats='mean').values]], 
             fmt='o', color='purple', capsize=5,  label=f"IC {IC} [1995-2014] = {(m_diff_30_TM.sel(stats='upper').values - m_diff_30_TM.sel(stats='lower').values):.2f} ZJ")

ax.hlines(m_diff_30_DL.sel(stats='mean').values, 1995, 2014, colors='lightgreen', linewidth=3)
ax.errorbar(x=2005, y=m_diff_30_DL.sel(stats='mean').values , 
             yerr=[[m_diff_30_DL.sel(stats='mean').values  - m_diff_30_DL.sel(stats='lower').values ], [m_diff_30_DL.sel(stats='upper').values  - m_diff_30_DL.sel(stats='mean').values]], 
             fmt='o', color='lightgreen', capsize=5, label=f"IC {IC} [1995-2014] = {(m_diff_30_DL.sel(stats='upper').values - m_diff_30_DL.sel(stats='lower').values):.2f} ZJ")

ax.axhline(0, color='green', linewidth=3)

ax.set_xlabel("Time (y)", fontsize=18)
ax.set_ylabel(f"Anomaly {var} / {ref_start}_{ref_end} ({unit})", fontsize = 18)
ax.set_title("[Hist_dd+3000 - Hist_dd+1000] depending on dedrift method", fontsize=22 , fontweight = "bold")
ax.tick_params(axis='both', labelsize=15)
ax.axvspan(ref_start, ref_end, color='grey', alpha=0.2)
ax.axvspan(p_start, p_end, color='grey', alpha=0.2)
ax.legend(loc='lower center', bbox_to_anchor=(0.3, 0.7), ncol=1, fontsize=15)
ax.grid()

plt.savefig(FIG_DIR / "figure_10.pdf", bbox_inches="tight")
plt.savefig(FIG_DIR / "figure_10.png", bbox_inches="tight")
