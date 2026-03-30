import xarray as xr
from dask.diagnostics import ProgressBar
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
import datetime 


plt.rc('font', size=12.5) #controls default text size
plt.rcParams['figure.dpi']=300 #controls fig resolution
edir = "/cnrm/cmip/CMIP6/CMIP/CNRM-CERFACS/CNRM-CM6-1/piControl/r1i1p1f2/Amon/tas/gr/latest"
edir_area_oce= '/cnrm/cmip/CMIP6/CMIP/CNRM-CERFACS/CNRM-CM6-1/piControl/r1i1p1f2/Ofx/areacello/gn/latest' 
outdir = "/cnrm/ioga/Users/seguy/hist_ensbl/pic_global"

SAT = xr.open_mfdataset(f"{outdir}/tas_gb_pic_3000y.nc", use_cftime = True)
vdata_sat = SAT.__xarray_dataarray_variable__

SST = xr.open_mfdataset(f"{outdir}/tos_gb_pic_3000y.nc", use_cftime = True)
vdata_sst = SST.__xarray_dataarray_variable__

MOC = xr.open_mfdataset(f"{outdir}/moc_pic_gb_3000y.nc", use_cftime = True)
vdata_moc = MOC

DOT = xr.open_mfdataset(f"{outdir}/dot_pic_gb_3000y.nc", use_cftime = True)
vdata_dot = DOT.thetao

TOA = xr.open_mfdataset(f"{outdir}/toa_sum_pic_gb_3000y.nc", use_cftime = True)
vdata_toa = TOA.__xarray_dataarray_variable__
vdata_toa['year'] = np.linspace(0,3000,3001)


slopes, years = rolling_trend_np(vdata_sat, window=1500)
slopes_100y_sat = slopes * 100

slopes, years = rolling_trend_np(vdata_sst, window=1500)
slopes_100y_sst = slopes * 100

mean_100y_toa = rolling_mean_np(vdata_toa.values, window=1500)
slopes, years = rolling_trend_np(vdata_toa.values, window=1500)
slopes_100y_toa = slopes * 100

slopes, years = rolling_trend_np(vdata_dot, window=1500)
slopes_100y_dot = slopes * 100

slopes, years = rolling_trend_np(vdata_moc.AMOC, window=1500)
slopes_1000y_amoc = slopes * 1000
slopes, years = rolling_trend_np(vdata_moc.PMOC, window=1500)
slopes_1000y_pmoc = slopes * 1000
slopes, years = rolling_trend_np(vdata_moc.SOMOC, window=1500)
slopes_1000y_somoc = slopes * 1000
slopes, years = rolling_trend_np(vdata_moc.GMOC, window=1500)
slopes_1000y_gmoc = slopes * 1000


fig, axs = plt.subplots(3, 1, figsize=(15, 15))

panel_labels = ["a)", "b)", "c)"]

# ============================================================
# Premier sous-graphe
# ============================================================
axs[0].plot(np.linspace(1500,2999,1500), vdata_toa.isel(year=slice(1500,3000)), color='grey')
axs[0].plot(np.linspace(1499,3001,1502), mean_100y_toa, color='black', label='TOA last 1000y mean')
axs[0].hlines(0.1, 1500, 3000, color="red", label="< 0.1 W m⁻²", linewidth=5)
axs[0].hlines(0.1, 1500, 3000, color="green", label="< 0.1 W m⁻²", linewidth=2)
axs[0].set_ylabel("TOA flux (W m⁻²)", fontsize=20, fontweight="bold")
axs[0].tick_params(axis='both', which='major', labelsize=20)
axs[0].set_title("Radiative balance", fontsize=25, fontweight="bold")
axs[0].grid(True)
axs[0].legend(loc='lower right', fontsize=15)

# --- label a) (hors cadre)
axs[0].text(
    -0.08, 1.02, panel_labels[0],
    transform=axs[0].transAxes,
    ha="right", va="bottom",
    fontsize=22, fontweight="bold",
    clip_on=False
)

# ============================================================
# Deuxième sous-graphe
# ============================================================
axs[1].plot(years, slopes_100y_sat, color='black', label='SAT')
axs[1].plot(years, slopes_100y_sst, color='grey', label='SST')
axs[1].plot(years, slopes_100y_dot, '--', color='grey', label='DOT')
axs[1].hlines(0.05,1500, 3000, color="red", label="< 0.05 °C century⁻¹", linewidth=5)
axs[1].hlines(0.1, 1500, 3000, color="green", label="< 0.1 °C century⁻¹", linewidth=2)
axs[1].set_ylabel("1000y Trends (°C century⁻¹)", fontsize=20, fontweight="bold")
axs[1].tick_params(axis='both', which='major', labelsize=20)
axs[1].set_title("Temperature", fontsize=25, fontweight="bold")
axs[1].grid(True)
axs[1].legend(fontsize=15)

# --- label b) (hors cadre)
axs[1].text(
    -0.08, 1.02, panel_labels[1],
    transform=axs[1].transAxes,
    ha="right", va="bottom",
    fontsize=22, fontweight="bold",
    clip_on=False
)

# ============================================================
# Troisième sous-graphe
# ============================================================
axs[2].plot(years, slopes_1000y_amoc, color='black', label='AMOC')
axs[2].plot(years, slopes_1000y_pmoc, '--', color='black', label='PMOC')
axs[2].plot(years, slopes_1000y_gmoc, color='black', linewidth=3, label='GMOC')
axs[2].plot(years, slopes_1000y_somoc, color='grey', label='SOMOC')
axs[2].hlines(1, 1500, 3000, color="green", linewidth=2)
axs[2].hlines(-1, 1500, 3000, color="green", linewidth=2)
axs[2].fill_between(years, 1, -1, color='lightgreen', alpha=0.5, label='< 1 Sv millenial⁻¹')
axs[2].set_xlabel("Simulated years", fontsize=20)
axs[2].tick_params(axis='both', which='major', labelsize=20)
axs[2].set_ylabel("1000y Trends (Sv millenial⁻¹)", fontsize=20, fontweight="bold")
axs[2].set_title("MOC", fontsize=25, fontweight="bold")
axs[2].grid(True)
axs[2].legend(fontsize=15)

# --- label c) (hors cadre)
axs[2].text(
    -0.08, 1.02, panel_labels[2],
    transform=axs[2].transAxes,
    ha="right", va="bottom",
    fontsize=22, fontweight="bold",
    clip_on=False
)

plt.tight_layout()













