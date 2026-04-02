
"""
Figure 10: AMOC and OHC changes correlation in atlantic subpolar gyre
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
INTERMEDIATE_DIR = PROJECT_ROOT / "data/data_plot"
METADATA_DIR =  PROJECT_ROOT / "data/metadata"


boot_3000_ohc = xr.open_mfdataset(INTERMEDIATE_DIR / "boot_3000_ohc.nc").__xarray_dataarray_variable__
boot_3000_amoc = xr.open_dataset(INTERMEDIATE_DIR / "boot_3000_amoc.nc").__xarray_dataarray_variable__
boot_1000_ohc = xr.open_mfdataset(INTERMEDIATE_DIR / "boot_1000_ohc.nc").__xarray_dataarray_variable__
boot_1000_amoc = xr.open_dataset(INTERMEDIATE_DIR / "boot_1000_amoc.nc").__xarray_dataarray_variable__

AMOC_cor_trends_1000 = xr.open_dataset(INTERMEDIATE_DIR / "AMOC_cor_trends_1000.nc").amoc
AMOC_cor_trends_3000 = xr.open_dataset(INTERMEDIATE_DIR / "AMOC_cor_trends_3000.nc").amoc
ohc_2D_cor_trends_gsp_1000 = xr.open_dataset(INTERMEDIATE_DIR / "ohc_2D_cor_trends_gsp_1000.nc").__xarray_dataarray_variable__
ohc_2D_cor_trends_gsp_3000 = xr.open_dataset(INTERMEDIATE_DIR / "ohc_2D_cor_trends_gsp_3000.nc").__xarray_dataarray_variable__

ohc_2D_cor_trends_1000 = xr.open_dataset(INTERMEDIATE_DIR / "ohc_2D_cor_trends_1000.nc").__xarray_dataarray_variable__
ohc_2D_cor_trends_3000 = xr.open_dataset(INTERMEDIATE_DIR / "ohc_2D_cor_trends_3000.nc").__xarray_dataarray_variable__


gsubpolmsk = xr.open_mfdataset(METADATA_DIR / "gsubpolmsk.nc").gsubpolmsk
ohc_2D_cor_trends_diff = ((ohc_2D_cor_trends_3000.mean(dim='ensemble') - ohc_2D_cor_trends_1000.mean(dim='ensemble'))*(1e12))*gsubpolmsk


# MEAN
ohc_cor_mean_1000 = boot_1000_ohc.sel(stats = "mean")
amoc_cor_mean_1000 = boot_1000_amoc.sel(stats = "mean")

ohc_cor_mean_3000 = boot_3000_ohc.sel(stats = "mean")
amoc_cor_mean_3000 = boot_3000_amoc.sel(stats = "mean")


# IC 90% 
ohc_cor_inf_1000 = boot_1000_ohc.sel(stats = "lower")
amoc_cor_inf_1000 = boot_1000_amoc.sel(stats = "lower")

ohc_cor_inf_3000 = boot_3000_ohc.sel(stats = "lower")
amoc_cor_inf_3000 = boot_3000_amoc.sel(stats = "lower")

ohc_cor_sup_1000 = boot_1000_ohc.sel(stats = "upper")
amoc_cor_sup_1000 = boot_1000_amoc.sel(stats = "upper")

ohc_cor_sup_3000 = boot_3000_ohc.sel(stats = "upper")
amoc_cor_sup_3000 = boot_3000_amoc.sel(stats = "upper")


#CORRELATION 

x_1000 = np.array([val[()] for val in AMOC_cor_trends_1000])
x_3000 = np.array([val[()] for val in AMOC_cor_trends_3000])

y_1000 = ohc_2D_cor_trends_gsp_1000
y_3000 = ohc_2D_cor_trends_gsp_3000

correlation_coefficient_1000, _ = stats.pearsonr(x_1000, y_1000)
slope_1000, intercept_1000, r_value_1000, p_value_1000, std_err_1000 = stats.linregress(x_1000, y_1000)
line_1000 = slope_1000 * x_1000 + intercept_1000

correlation_coefficient_3000, _ = stats.pearsonr(x_3000, y_3000)
slope_3000, intercept_3000, r_value_3000, p_value_3000, std_err_3000 = stats.linregress(x_3000, y_3000)
line_3000 = slope_3000 * x_3000 + intercept_3000


#TRACE

# erreurs 3000
xerr_3000 = [[amoc_cor_mean_3000 - amoc_cor_inf_3000],
             [amoc_cor_sup_3000 - amoc_cor_mean_3000]]

yerr_3000 = [[ohc_cor_mean_3000 - ohc_cor_inf_3000],
             [ohc_cor_sup_3000 - ohc_cor_mean_3000]]

# erreurs 1000
xerr_1000 = [[amoc_cor_mean_1000 - amoc_cor_inf_1000],
             [amoc_cor_sup_1000 - amoc_cor_mean_1000]]

yerr_1000 = [[ohc_cor_mean_1000 - ohc_cor_inf_1000],
             [ohc_cor_sup_1000 - ohc_cor_mean_1000]]



# =========================
# Figure map
# =========================

fig = plt.figure(figsize=(5.2, 3.6), dpi=300)
ax = plt.axes(projection=ccrs.PlateCarree())

# Adapter les noms de coordonnées si besoin
lon_name = "lon" if "lon" in ohc_2D_cor_trends_diff.coords else "longitude"
lat_name = "lat" if "lat" in ohc_2D_cor_trends_diff.coords else "latitude"

lon = ohc_2D_cor_trends_diff[lon_name]
lat = ohc_2D_cor_trends_diff[lat_name]

# Colormap / limites
vmin = -7
vmax = 0

pcm = ax.pcolormesh(
    lon,
    lat,
    ohc_2D_cor_trends_diff,
    transform=ccrs.PlateCarree(),
    cmap="viridis",
    vmin=vmin,
    vmax=vmax,
    shading="auto"
)

# Côte / fond
ax.coastlines(linewidth=0.6)
ax.add_feature(cfeature.LAND, facecolor="white", edgecolor="none", zorder=10)

# Zoom proche de ce qu'on voit sur ton image
ax.set_extent([-85, 25, 20, 80], crs=ccrs.PlateCarree())

# Retirer le cadre/ticks si tu veux juste l'inset
ax.set_xticks([])
ax.set_yticks([])
ax.set_frame_on(False)

# Colorbar
cbar = plt.colorbar(pcm, ax=ax, fraction=0.05, pad=0.03)
cbar.set_label(
    "OHC gain diff mask (GJ m$^{-2}$)",
    rotation=270,
    labelpad=15,
    fontsize=10
)
cbar.ax.tick_params(labelsize=8)
cbar.ax.tick_params(labelsize=8)


plt.savefig(FIG_DIR / "figure_10_2.pdf", bbox_inches="tight", pad_inches=0.02)
plt.savefig(FIG_DIR / "figure_10_2.png", bbox_inches="tight", dpi=300)





# =========================
# Figure correlation
# =========================


plt.figure(figsize=(15, 10))
plt.scatter(x_3000, y_3000, color='blue', label='hist+3000')
plt.plot(x_3000, line_3000, color='blue', label='Linear regression hist_dd+3000', linewidth=2)
plt.plot(amoc_cor_mean_3000, ohc_cor_mean_3000, 'x', color = 'blue', label = 'Ensemble mean', markersize = 30)
plt.scatter(x_1000, y_1000, color='orange', label='hist+1000')
plt.plot(x_1000, line_1000, color='orange', label='Linear regression hist_dd+1000', linewidth=2)
plt.plot(amoc_cor_mean_1000, ohc_cor_mean_1000, 'x', color = 'orange', label = 'Ensemble mean', markersize = 30)
plt.errorbar(amoc_cor_mean_3000, ohc_cor_mean_3000,
             xerr=xerr_3000, yerr=yerr_3000,
             fmt='x', color='blue', markersize=15,
             capsize=5)

plt.errorbar(amoc_cor_mean_1000, ohc_cor_mean_1000,
             xerr=xerr_1000, yerr=yerr_1000,
             fmt='x', color='orange', markersize=15,
             capsize=5)
plt.text(min(x_1000)+4, min(y_1000.values),
         f'r = {r_value_1000:.2f}, p = {p_value_1000:.2e}',
         fontsize=15, color='orange', fontweight='bold', verticalalignment='top')

plt.text(min(x_3000)+2, min(y_3000.values),
         f'r = {r_value_3000:.2f}, p = {p_value_3000:.2e}',
         fontsize=15, color='blue', fontweight='bold', verticalalignment='top')

plt.xlabel('AMOC historical dedrifted change (Sv)', fontsize=20)
plt.ylabel('OHC historical dedrifted gains (ZJ)', fontsize=20)
plt.legend()
plt.grid()


plt.savefig(FIG_DIR / "figure_10_1.pdf", bbox_inches="tight")
plt.savefig(FIG_DIR / "figure_10_1.png", bbox_inches="tight", dpi=300)


# ============================================================
# FINAL FIGURE
# ============================================================

main_fig_path = FIG_DIR / "figure_10_1.png"  
inset_fig_path = FIG_DIR / "figure_10_2.png"     
output_path = FIG_DIR / "figure_10.pdf"


main_img = Image.open(main_fig_path).convert("RGBA")
inset_img = Image.open(inset_fig_path).convert("RGBA")
inset_pos = [0.68, 0.08, 0.32, 0.28]

w, h = main_img.size
dpi = 100
fig = plt.figure(figsize=(w / dpi, h / dpi), dpi=dpi)

ax_main = fig.add_axes([0, 0, 1, 1])
ax_main.imshow(main_img)
ax_main.axis("off")

# Axe inset
ax_inset = fig.add_axes(inset_pos)
ax_inset.imshow(inset_img)
ax_inset.axis("off")


plt.savefig(FIG_DIR / "figure_10.pdf", bbox_inches="tight")
plt.savefig(FIG_DIR / "figure_10.png", bbox_inches="tight")