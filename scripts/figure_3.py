
import numpy as np
import xarray as xr
import statsmodels.api as sm
from scipy.stats import norm
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors
from dask.diagnostics import ProgressBar
import datetime 
from scipy import stats
from matplotlib.cm import get_cmap
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.patches as mpatches


edir_pic = "/cnrm/ioga/Users/seguy/hist_OHC"
edir_area_oce=  "/home/seguyr/Documents/TT/areacello" #'/cnrm/cmip/CMIP6/CMIP/CNRM-CERFACS/CNRM-CM6-1/piControl/r1i1p1f2/Ofx/areacello/gn/latest' 


var = 'ohc'
echelle = 1e12
Cp = 3991.867 # J/(kg*K)
rho = 1026 # kg/m³
unit = 'GJ/m2'
layer_name = '0_btm'
n_boot = 1000 
span_width = 25
qt_inf = 0.05
qt_sup = 0.95
IC = '90'


# === Coordonnées grille nemo ===
farea_oce = xr.open_mfdataset("%s/*.nc"%edir_area_oce) 
area_oce = farea_oce['areacello']
area_tot = 3.626975e+14
area_poids = area_oce/area_tot

def correc_nemo_lon(lon) :
    after_discont = ~(lon.diff("x", label="upper") > 0).cumprod("x").astype(bool)
    lon[:,1:] = lon[:,1:] + after_discont*360.
    return lon
lon2d_oce = correc_nemo_lon(area_oce.lon)
lat2d_oce = area_oce.lat


ohc_pic_gb = xr.open_mfdataset("%s/2D_ohc_pic_global_0_btm.nc" % (edir_pic)).__xarray_dataarray_variable__ * 1e12
ohc_pic_gb['time'] = np.linspace(0,3000,3000)
ohc_1000 = ohc_pic_gb.isel(time=slice(0,1000))
ohc_400 =  ohc_pic_gb.isel(time=slice(2600,3000))



fit_map_1000 = fit_map_hac(ohc_1000, conf=90)
slope_1000 = fit_map_1000["slope"]*100
ci_low_1000 = fit_map_1000["ci_low"]*100
ci_high_1000 = fit_map_1000["ci_high"]*100
pvalue_1000 = fit_map_1000["pvalue"]*100


fit_map_400 = fit_map_hac(ohc_400, conf=90)
slope_400 = fit_map_400["slope"]*100
ci_low_400 = fit_map_400["ci_low"]*100
ci_high_400 = fit_map_400["ci_high"]*100
pvalue_400 = fit_map_400["pvalue"]*100


# -----------------------------------------------------------------------------
# Visualisation
# -----------------------------------------------------------------------------

VMIN, VMAX = -3, 3
WHITE_WIDTH = 0.01 * (VMAX - VMIN)

bounds = np.linspace(VMIN, VMAX, 256)
colors = plt.cm.RdBu_r(np.linspace(0, 1, 256))
i0_low = np.argmin(np.abs(bounds + WHITE_WIDTH))
i0_high = np.argmin(np.abs(bounds - WHITE_WIDTH))
colors[i0_low:i0_high] = [1, 1, 1, 1]

CMAP_CUSTOM = mcolors.ListedColormap(colors)
NORM = mcolors.TwoSlopeNorm(vmin=VMIN, vcenter=0, vmax=VMAX)


# -----------------------------
# FIGURE 
# -----------------------------
proj = ccrs.Robinson(central_longitude=0)

fig = plt.figure(figsize=(15, 5))
gs = fig.add_gridspec(
    nrows=2, ncols=2,
    height_ratios=[1, 0.05],
    hspace=0.1, wspace=0.03
)

axes = [
    fig.add_subplot(gs[0, 0], projection=proj),  # gauche
    fig.add_subplot(gs[0, 1], projection=proj)   # droite
]

cax = fig.add_subplot(gs[1, :])  # colorbar sur toute la largeur

cf = plot_panel(axes[0], "Pic+1000", "a)", CMAP_CUSTOM, NORM)
plot_panel(axes[1], "Pic+3000", "b)", CMAP_CUSTOM, NORM)

cbar = fig.colorbar(
    cf, cax=cax,
    orientation='horizontal',
    extend='both',
    extendfrac='auto'
)

cbar.ax.xaxis.set_label_position('bottom')
cbar.ax.xaxis.set_ticks_position('bottom')
cbar.set_label("OHC trend (GJ m$^{-2}$century$^{-1}$)", fontsize=20, labelpad=20)
cbar.ax.tick_params(labelsize=20, length=7, width=1.2)
