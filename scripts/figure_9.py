"""
Figure 9: obs VS model comparison for SAT LAND and SST
====================================================

"""

from pathlib import Path
import sys
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import matplotlib.gridspec as gridspec


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
    COLORS, 
    plot_row
)


t_m = np.linspace(1850,2014,165)
t_o =  np.linspace(1940,2014,75)


scalars_ds = xr.open_dataset(INTERMEDIATE_DIR/"scalars.nc")

m_boot_1000 = xr.open_dataset(INTERMEDIATE_DIR/"m_boot_1000.nc").__xarray_dataarray_variable__
m_boot_3000 = xr.open_dataset(INTERMEDIATE_DIR/"m_boot_3000.nc").__xarray_dataarray_variable__
boot_1000 = xr.open_dataset(INTERMEDIATE_DIR/"boot_1000.nc").__xarray_dataarray_variable__
boot_3000 = xr.open_dataset(INTERMEDIATE_DIR/"boot_3000.nc").__xarray_dataarray_variable__

m_boot_1000_2 = xr.open_dataset(INTERMEDIATE_DIR/"m_boot_1000_2.nc").__xarray_dataarray_variable__
m_boot_3000_2 = xr.open_dataset(INTERMEDIATE_DIR/"m_boot_3000_2.nc").__xarray_dataarray_variable__
boot_1000_2 = xr.open_dataset(INTERMEDIATE_DIR/"boot_1000_2.nc").__xarray_dataarray_variable__
boot_3000_2 = xr.open_dataset(INTERMEDIATE_DIR/"boot_3000_2.nc").__xarray_dataarray_variable__

obs_mean = xr.open_dataset(INTERMEDIATE_DIR/"obs_mean.nc").sst
obs_ic = xr.open_dataset(INTERMEDIATE_DIR/"obs_ic.nc").sst
obs_mean_2 = xr.open_dataset(INTERMEDIATE_DIR/"obs_mean_2.nc").tasland
obs_ic_2 = xr.open_dataset(INTERMEDIATE_DIR/"obs_ic_2.nc").tasland


error_tot = scalars_ds["error_tot"].item()
error_tot_2 = scalars_ds["error_tot_2"].item()
tas_gain = scalars_ds["tas_gain"].item()
tas_gain_2 = scalars_ds["tas_gain_2"].item()

gain_tot = m_boot_1000.sel(stats='mean')
qsup_tot = m_boot_1000.sel(stats='upper')
qinf_tot = m_boot_1000.sel(stats='lower')

gain_3000 = m_boot_3000.sel(stats='mean')
qsup_3000 = m_boot_3000.sel(stats='upper')
qinf_3000 = m_boot_3000.sel(stats='lower')


gain_tot_2 = m_boot_1000_2.sel(stats='mean')
qsup_tot_2 = m_boot_1000_2.sel(stats='upper')
qinf_tot_2 = m_boot_1000_2.sel(stats='lower')

gain_3000_2 = m_boot_3000_2.sel(stats='mean')
qsup_3000_2 = m_boot_3000_2.sel(stats='upper')
qinf_3000_2 = m_boot_3000_2.sel(stats='lower')



# ============================================================
# FIGURE 2x2
# ============================================================
fig = plt.figure(figsize=(35, 25))
gs = gridspec.GridSpec(
    2, 2, figure=fig,
    width_ratios=[4, 2],
    height_ratios=[1, 1],
    wspace=0.02, 
    hspace=0.25
)

ax_a = fig.add_subplot(gs[0, 0])  # variable 1 - time series
ax_b = fig.add_subplot(gs[0, 1])  # variable 1 - bar plot
ax_c = fig.add_subplot(gs[1, 0])  # variable 2 - time series
ax_d = fig.add_subplot(gs[1, 1])  # variable 2 - bar plot


# ============================================================
# LIGNE 1 : VARIABLE 1
# ============================================================
plot_row(
    ax_ts=ax_a,
    ax_bar=ax_b,
    t_m=t_m,
    boot_1000=boot_1000,
    boot_3000=boot_3000,
    t_o=t_o,
    obs_mean=obs_mean,
    obs_ic=obs_ic,
    tas_gain=tas_gain,
    error_tot=error_tot,
    gain_tot=gain_tot,
    qinf_tot=qinf_tot,
    qsup_tot=qsup_tot,
    gain_3000=gain_3000,
    qinf_3000=qinf_3000,
    qsup_3000=qsup_3000,
    ylabel_ts="Sea Surface Temperature (°C)",
    ylabel_bar="Change over 1958–2014 (°C)",
    label_obs = "ERA5c, AMIP-1-1-10, AMIP-ERSST5-1-0, AMIP-Had1p1-1-0",
    legend_y=-0.1
)

# ============================================================
# LIGNE 2 : VARIABLE 2
# Remplacer ici par tes variables 2
# ============================================================
plot_row(
    ax_ts=ax_c,
    ax_bar=ax_d,
    t_m=t_m,
    boot_1000=boot_1000_2,
    boot_3000=boot_3000_2,
    t_o=t_o,
    obs_mean=obs_mean_2,
    obs_ic=obs_ic_2,
    tas_gain=tas_gain_2,
    error_tot=error_tot_2,
    gain_tot=gain_tot_2,
    qinf_tot=qinf_tot_2,
    qsup_tot=qsup_tot_2,
    gain_3000=gain_3000_2,
    qinf_3000=qinf_3000_2,
    qsup_3000=qsup_3000_2,
    ylabel_ts="Land 2-m Temperature (°C)",
    ylabel_bar="Change over 1958–2014 (°C)",
    label_obs = "ERA5c, BEST, CRU-TS4.03",
    legend_y=-0.1
)

fig.text(
    0.05, 0.90, "a)",
    ha="left", va="top",
    fontsize=40, fontweight="bold"
)

fig.text(
    0.05, 0.48, "b)",
    ha="left", va="top",
    fontsize=40, fontweight="bold"
)


plt.savefig(FIG_DIR / "figure_9.pdf", bbox_inches="tight")
plt.savefig(FIG_DIR / "figure_9.png", bbox_inches="tight")






