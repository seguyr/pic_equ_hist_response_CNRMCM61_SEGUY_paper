"""
Figure 12: Dedrift and sampling results robustess
"""

from pathlib import Path
import sys
import matplotlib.pyplot as plt
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



ref_start = 1850
ref_end = 1900
p_start = 1995
p_end = 2014
x_gb, y_gb = 0.7, 0.84
IC = "90%"
var = 'ohc'
unit = 'ZJ'


diff_hist_cor_tot = xr.open_mfdataset(INTERMEDIATE_DIR / "diff_hist_cor_tot.nc").__xarray_dataarray_variable__
m_diff_hist_cor_tot = xr.open_mfdataset(INTERMEDIATE_DIR / "m_diff_hist_cor_tot.nc").__xarray_dataarray_variable__
diff_hist_cor_tot_TM = xr.open_mfdataset(INTERMEDIATE_DIR / "diff_hist_cor_tot_TM.nc").__xarray_dataarray_variable__
m_diff_hist_cor_tot_TM = xr.open_mfdataset(INTERMEDIATE_DIR / "m_diff_hist_cor_tot_TM.nc").__xarray_dataarray_variable__
diff_hist_cor_1000_3 = xr.open_mfdataset(INTERMEDIATE_DIR / "diff_hist_cor_1000_3.nc").__xarray_dataarray_variable__
m_diff_hist_cor_1000_3 = xr.open_mfdataset(INTERMEDIATE_DIR / "m_diff_hist_cor_1000_3.nc").__xarray_dataarray_variable__




# === Tracé ===
y = y_gb
x = x_gb
time = diff_hist_cor_1000_3.time + ref_start
fig, ax = plt.subplots(1, 1, figsize=(15, 10), sharey= False)


ax.plot(time, diff_hist_cor_tot_TM.sel(stats='mean').values, color="purple", label='Time Matching dedrift with 29 members')
ax.fill_between(time, diff_hist_cor_tot_TM.sel(stats='lower').values, diff_hist_cor_tot_TM.sel(stats='upper').values, color="purple", alpha=0.7)

ax.hlines(m_diff_hist_cor_tot_TM.sel(stats='mean').values, 1995, 2014, colors='purple',  linewidth=3)
ax.errorbar(x=2005, y=m_diff_hist_cor_tot_TM.sel(stats='mean').values, 
             yerr=[[m_diff_hist_cor_tot_TM.sel(stats='mean').values  - m_diff_hist_cor_tot_TM.sel(stats='lower').values ], [m_diff_hist_cor_tot_TM.sel(stats='upper').values  - m_diff_hist_cor_tot_TM.sel(stats='mean').values]], 
             fmt='o',linewidth = 7, color='purple', capsize=15, label=f"IC {IC} [1995-2014] = {(m_diff_hist_cor_tot_TM.sel(stats='upper').values - m_diff_hist_cor_tot_TM.sel(stats='lower').values):.2f} ZJ")


ax.plot(time, diff_hist_cor_tot.sel(stats='mean').values, color="red", label='Linear Regression dedrift with 29 members')
ax.fill_between(time, diff_hist_cor_tot.sel(stats='lower').values, diff_hist_cor_tot.sel(stats='upper').values, color="red", alpha=0.4)

ax.hlines(m_diff_hist_cor_tot.sel(stats='mean').values, 1995, 2014, colors='red', linewidth=3)
ax.errorbar(x=2005, y=m_diff_hist_cor_tot.sel(stats='mean').values, 
             yerr=[[m_diff_hist_cor_tot.sel(stats='mean').values  - m_diff_hist_cor_tot.sel(stats='lower').values ], [m_diff_hist_cor_tot.sel(stats='upper').values  - m_diff_hist_cor_tot.sel(stats='mean').values]], 
             fmt='o', linewidth = 5, color='red', capsize=15, label=f"IC {IC} [1995-2014] = {(m_diff_hist_cor_tot.sel(stats='upper').values - m_diff_hist_cor_tot.sel(stats='lower').values):.2f} ZJ")

ax.plot(time, diff_hist_cor_1000_3.sel(stats='mean').values, color="yellow", label='Linear Regression dedrift with 10 members')
ax.fill_between(time, diff_hist_cor_1000_3.sel(stats='lower').values, diff_hist_cor_1000_3.sel(stats='upper').values, color="yellow", alpha=0.6)

ax.hlines(m_diff_hist_cor_1000_3.sel(stats='mean').values, 1995, 2014, colors='yellow', linewidth=3)
ax.errorbar(x=2005, y=m_diff_hist_cor_1000_3.sel(stats='mean').values, 
             yerr=[[m_diff_hist_cor_1000_3.sel(stats='mean').values  - m_diff_hist_cor_1000_3.sel(stats='lower').values ], [m_diff_hist_cor_1000_3.sel(stats='upper').values  - m_diff_hist_cor_1000_3.sel(stats='mean').values]], 
             fmt='o', linewidth = 3, color='yellow', capsize=15,  label=f"IC {IC} [1995-2014] = {(m_diff_hist_cor_1000_3.sel(stats='upper').values - m_diff_hist_cor_1000_3.sel(stats='lower').values):.2f} ZJ")

ax.axhline(0, color='green', linewidth=3)

ax.set_xlabel("Time (y)", fontsize=18)
ax.set_ylabel(f"Anomaly {var} / {ref_start}_{ref_end} ({unit})", fontsize = 18)
ax.set_title("[Hist_dd+3000 - Hist_dd+1000]", fontsize=22 , fontweight = "bold")
ax.tick_params(axis='both', labelsize=15)
ax.axvspan(ref_start, ref_end, color='grey', alpha=0.2)
ax.axvspan(p_start, p_end, color='grey', alpha=0.2)
ax.legend(loc='lower center', bbox_to_anchor=(0.3, 0.7), ncol=1, fontsize=15)
ax.grid()
plt.savefig(FIG_DIR / "figure_12.pdf", bbox_inches="tight")
plt.savefig(FIG_DIR / "figure_12.png", bbox_inches="tight")