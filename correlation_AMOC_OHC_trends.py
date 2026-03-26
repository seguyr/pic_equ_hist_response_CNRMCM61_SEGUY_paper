#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  2 15:44:03 2025

@author: seguyr
"""
#calculer l'amoc de chaque hist du LR

import xarray as xr
from dask.diagnostics import ProgressBar
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.colors as mcolors
from scipy import stats


plt.rc('font', size=12.5) #controls default text size
plt.rcParams['figure.dpi']=300 #controls fig resolution
edir_hist = "/cnrm/cmip/CMIP6/CMIP/CNRM-CERFACS/CNRM-CM6-1/historical"
outdir = "/cnrm/ioga/Users/seguy/hist_ensbl" 
edir = "/cnrm/ioga/Users/seguy/hist_ens" 
edir_area= '/cnrm/cmip/CMIP6/CMIP/CNRM-CERFACS/CNRM-CM6-1/piControl/r1i1p1f2/fx/areacella/gr/latest' 
edir_area_oce=  "/home/seguyr/Documents/TT/areacello" # '/cnrm/cmip/CMIP6/CMIP/CNRM-CERFACS/CNRM-CM6-1/piControl/r1i1p1f2/Ofx/areacello/gn/latest' 
edir_pic = "/cnrm/cmip/CMIP6/CMIP/CNRM-CERFACS/CNRM-CM6-1/picontrol"
outdir_TT = "Documents/TT/hist_ensbl/fig_papier"


#tracé
print('selection areacello')
farea_oce = xr.open_mfdataset("%s/*.nc"%edir_area_oce) 
area_oce = farea_oce['areacello']


echelle_amoc = (10**6)  #sverdrup
rho = 1026 # kg/m³ 
echelle_ohc = (10**9)  #GJ
area_tot = 3.626975e+14  #ocean
Cp = 3991.867 # J/(kg*K)
rho = 1026 # kg/m³ 
var = 'ohc'
n_boot = 1000
qt_inf = 0.05
qt_sup = 0.95

liste_hist = [12,13,14,15,16,17,18,28,29,30,31,32,33,34,35,36,37,38,39,40]
liste_autre = [2,3,4,5,6,7,8,9,10,11,19,20,21,22,23,24,25,26,27]
liste_1000 = [12,13,14,15,16,17,18,28,29,30]
liste_tot = [k for k in range(2,31)]



def bootstrap(arr, n_boot=n_boot, seed=0):
    if np.all(np.isnan(arr)):
        return np.array([np.nan, np.nan, np.nan, np.nan, np.nan], dtype=np.float32)

    # Générateur pseudo-aléatoire
    rng = np.random.default_rng(seed)

    n_ens = arr.sizes["ensemble"]
    ens_values = arr.ensemble.values

    print("sample")
    sampled_ens = rng.choice(
        ens_values,
        size=(n_boot, n_ens),
        replace=True
    )  # (n_boot, n_ens)

    sampled_da = xr.DataArray(sampled_ens, dims=("boot", "ensemble"))

    print("resample")
    resampled = arr.sel(ensemble=sampled_da)

    print("mean")
    boot_means = resampled.mean(dim="ensemble")
    boot_means = boot_means.chunk(dict(boot=-1))
    
    print("stats")
    q_inf = boot_means.quantile(qt_inf, dim="boot")
    q_sup = boot_means.quantile(qt_sup, dim="boot")
    mean  = boot_means.mean(dim="boot")
    sigma = boot_means.std(dim="boot")

    p_left  = (boot_means <= 0).mean(dim="boot")
    p_right = (boot_means >= 0).mean(dim="boot")
    p_value = 2 * np.minimum(p_left, p_right)

    return np.array([q_inf, mean, q_sup, sigma, p_value], dtype=np.float32)


#AMOC index

amoc = []
amoc_pic = []
for k in liste_1000:
    print(k)
    ds = xr.open_mfdataset("%s/hist_%s/amoc/amoc_hist_r%s_165y_45N.nc"%(outdir,k, k))/(rho*echelle_amoc)
    amoc.append(ds)
    ds = xr.open_mfdataset("%s/hist_%s/amoc/amoc_pic_r%s_165y_45N.nc"%(outdir,k, k))/(rho*echelle_amoc)
    amoc_pic.append(ds)
AMOC_1000 = xr.concat(amoc, dim='ensemble')  
AMOC_pic_1000 = xr.concat(amoc_pic, dim='ensemble')  



amoc = []
amoc_pic = []
for k in liste_tot:
    print(k)
    ds = xr.open_mfdataset("%s/hist_%s/amoc/amoc_hist_r%s_165y_45N.nc"%(outdir,k, k))/(rho*echelle_amoc)
    amoc.append(ds)
    ds = xr.open_mfdataset("%s/hist_%s/amoc/amoc_pic_r%s_165y_45N.nc"%(outdir,k, k))/(rho*echelle_amoc)
    amoc_pic.append(ds)
AMOC_tot = xr.concat(amoc, dim='ensemble')  
AMOC_pic_tot = xr.concat(amoc_pic, dim='ensemble')  


amoc = []
amoc_pic = []
for k in range(31,41):
    print(k)
    ds = xr.open_mfdataset("%s/hist_%s/amoc/amoc_hist_r%s_165y_45N.nc"%(outdir,k, k))/(rho*echelle_amoc)
    amoc.append(ds)
    ds = xr.open_mfdataset("%s/hist_%s/amoc/amoc_pic_r%s_165y_45N.nc"%(outdir,k, k))/(rho*echelle_amoc)
    amoc_pic.append(ds)
    
AMOC_3000 = xr.concat(amoc, dim='ensemble')  
AMOC_pic_3000 = xr.concat(amoc_pic, dim='ensemble')  


edir_ohc = "/home/seguyr/Documents/TT/hist_ensbl/fig_papier"
output_path2 = "%s/pic_tot/amoc_pic_tot.nc"%(edir_ohc)
output_path3 = "%s/pic+3000/amoc_pic_3000.nc"%(edir_ohc)
output_path4 = "%s/hist_tot/amoc_hist_tot.nc"%(edir_ohc)
output_path5 = "%s/hist+3000/amoc_hist_3000.nc"%(edir_ohc)


with ProgressBar():
    AMOC_pic_tot.to_netcdf(path=output_path2, format="NETCDF4_CLASSIC")
    AMOC_pic_3000.to_netcdf(path=output_path3, format="NETCDF4_CLASSIC")
    AMOC_tot.to_netcdf(path=output_path4, format="NETCDF4_CLASSIC")
    AMOC_3000.to_netcdf(path=output_path5, format="NETCDF4_CLASSIC")



#OHC 
var = 'ohc'


OHC_2D_hist_tot = xr.open_mfdataset("/home/seguyr/%s/hist_tot/2D_ohc_hist_tot_0_btm.nc"%(outdir_TT),chunks={'time':120}, parallel=True, use_cftime=True)
OHC_2D_pic_tot = xr.open_mfdataset("/home/seguyr/%s/pic_tot/2D_ohc_pic_tot_0_btm.nc"%(outdir_TT),chunks={'time':120}, parallel=True, use_cftime=True)

OHC_2D_hist_3000 = xr.open_mfdataset("/home/seguyr/%s/hist+3000/2D_ohc_hist_3000_0_btm.nc"%(outdir_TT),chunks={'time':120}, parallel=True, use_cftime=True)
OHC_2D_pic_3000 = xr.open_mfdataset("/home/seguyr/%s/pic+3000/2D_ohc_pic_3000_0_btm.nc"%(outdir_TT),chunks={'time':120}, parallel=True, use_cftime=True)

OHC_2D_pic_tot['time'] = OHC_2D_hist_tot['time'] 

OHC_2D_pic_3000['time'] = OHC_2D_hist_3000['time'] 


#Ensembles corrigés
ohc_2D_cor_1000 = OHC_2D_hist_tot - OHC_2D_pic_tot
amoc_cor_1000 = AMOC_tot - AMOC_pic_tot

ohc_2D_cor_3000 = OHC_2D_hist_3000 - OHC_2D_pic_3000
amoc_cor_3000 = AMOC_3000 - AMOC_pic_3000


#tendances et moyennes
#ohc_2D_trends_1000 = (OHC_2D_1000.isel(time = slice(144,164)).mean(dim='time') - OHC_2D_1000.isel(time = slice(0,49)).mean(dim='time'))
ohc_2D_cor_trends_1000 = (ohc_2D_cor_1000.isel(time = slice(144,164)).mean(dim='time') - ohc_2D_cor_1000.isel(time = slice(0,49)).mean(dim='time'))
#AMOC_trends_1000 = (AMOC_1000.isel(year = slice(144,164)).mean(dim='year') - AMOC_1000.isel(year = slice(0,49)).mean(dim='year')).amoc.values
AMOC_cor_trends_1000 = (amoc_cor_1000.isel(year = slice(144,164)).mean(dim='year') - amoc_cor_1000.isel(year = slice(0,49)).mean(dim='year')).amoc
#AMOC_1900_1950_1000 = (AMOC_1000.sel(year=slice(1850,1900)).mean(dim='year')).amoc.values
#AMOC_cor_1900_1950_1000 = (amoc_cor_1000.sel(year=slice(1850,1900)).mean(dim='year')).amoc.values


# #tendances et moyennes
# ohc_2D_trends_autre = (OHC_2D_autre.isel(time = slice(144,164)).mean(dim='time') - OHC_2D_autre.isel(time = slice(0,49)).mean(dim='time'))
# ohc_2D_cor_trends_autre = (ohc_2D_cor_autre.isel(time = slice(144,164)).mean(dim='time') - ohc_2D_cor_autre.isel(time = slice(0,49)).mean(dim='time'))
# AMOC_trends_autre = (AMOC_autre.isel(year = slice(144,164)).mean(dim='year') - AMOC_autre.isel(year = slice(0,49)).mean(dim='year')).amoc.values
# AMOC_cor_trends_autre = (amoc_cor_autre.isel(year = slice(144,164)).mean(dim='year') - amoc_cor_autre.isel(year = slice(0,49)).mean(dim='year')).amoc.values
# AMOC_1900_1950_autre = (AMOC_autre.sel(year=slice(1850,1900)).mean(dim='year')).amoc.values
# AMOC_cor_1900_1950_autre = (amoc_cor_autre.sel(year=slice(1850,1900)).mean(dim='year')).amoc.values


#ohc_2D_trends_3000 = (OHC_2D_3000.isel(time = slice(144,164)).mean(dim='time') - OHC_2D_3000.isel(time = slice(0,49)).mean(dim='time'))
ohc_2D_cor_trends_3000 = (ohc_2D_cor_3000.isel(time = slice(144,164)).mean(dim='time') - ohc_2D_cor_3000.isel(time = slice(0,49)).mean(dim='time'))
#AMOC_trends_3000 = (AMOC_3000.isel(year = slice(144,164)).mean(dim='year') - AMOC_3000.isel(year = slice(0,49)).mean(dim='year')).amoc.values
AMOC_cor_trends_3000 = (amoc_cor_3000.isel(year = slice(144,164)).mean(dim='year') - amoc_cor_3000.isel(year = slice(0,49)).mean(dim='year')).amoc
#AMOC_1900_1950_3000 = (AMOC_3000.sel(year=slice(1850,1900)).mean(dim='year')).amoc.values
#AMOC_cor_1900_1950_3000 = (amoc_cor_3000.sel(year=slice(1850,1900)).mean(dim='year')).amoc.values



def bootstrap_2(arr1, arr2, n_boot, seed=0):
    if np.all(np.isnan(arr1)) or np.all(np.isnan(arr2)):
        return np.array([np.nan, np.nan, np.nan, np.nan, np.nan], dtype=np.float32)

    rng = np.random.default_rng(seed)

    n_ens1 = arr1.sizes["ensemble"]
    ens_values1 = arr1.ensemble.values
    
    n_ens2 = arr2.sizes["ensemble"]
    ens_values2 = arr2.ensemble.values
    
    sampled_ens_1 = rng.choice(ens_values1, size=(n_boot, n_ens1), replace=True)
    sampled_ens_2 = rng.choice(ens_values2, size=(n_boot, n_ens2), replace=True)  # indépendant du 1

    sampled_da_1 = xr.DataArray(sampled_ens_1, dims=("boot", "ensemble"))
    sampled_da_2 = xr.DataArray(sampled_ens_2, dims=("boot", "ensemble"))

    resampled_1 = arr1.sel(ensemble=sampled_da_1)
    resampled_2 = arr2.sel(ensemble=sampled_da_2)

    boot_means = resampled_2.mean(dim="ensemble") - resampled_1.mean(dim="ensemble")

    q_inf = boot_means.quantile(qt_inf, dim="boot")
    q_sup = boot_means.quantile(qt_sup, dim="boot")
    mean  = boot_means.mean(dim="boot")
    sigma = boot_means.std(dim="boot")

    p_left  = (boot_means <= 0).mean(dim="boot")
    p_right = (boot_means >= 0).mean(dim="boot")
    p_value = 2 * np.minimum(p_left, p_right)

    return np.array([q_inf, mean, q_sup, sigma, p_value], dtype=np.float32)


"""

bot_diff = bootstrap_2(ohc_2D_cor_trends_1000.__xarray_dataarray_variable__, ohc_2D_cor_trends_3000.__xarray_dataarray_variable__, n_boot)
bot_diff = xr.DataArray(bot_diff,dims=("stats", "x", "y"),coords={"stats": ["lower", "mean", "upper",  "sigma", "p_value"]})

mask = (bot_diff.where(bot_diff.sel(stats="p_value")<0.1).sel(stats="mean"))*(10**12)

ymin, ymax, xmin, xmax = 245, 265, 210, 235

valid = ~np.isnan(mask.values[xmin:xmax, ymin:ymax])
ix_local, iy_local = np.where(valid)

ix = ix_local + xmin
iy = iy_local + ymin

indices = np.column_stack([ix, iy])



# plot du mask
plt.figure()
ax = mask.plot()
quad = mask.plot()          
ax = quad.axes       
ax.scatter(iy, ix, s=10, c="red", alpha = 0.01)
plt.show()


# 1) prendre la variable du Dataset
da = ohc_2D_cor_trends_1000.__xarray_dataarray_variable__

# 2) construire un masque 2D sur (x,y)
tmpl = da.isel(ensemble=0)

arr = np.zeros(tmpl.shape, dtype=bool)
arr[ix, iy] = True

gsubpolmsk = xr.DataArray(arr, dims=tmpl.dims, coords=tmpl.coords)

# ton masque booléen
gsubpolmsk.name = "gsubpolmsk"

# option 1 : sauver directement le DataArray
gsubpolmsk.to_netcdf("Documents/TT/hist_ensbl/fig_papier/gsubpolmsk.nc")

"""

gsubpolmsk = xr.open_mfdataset("/home/seguyr/Documents/TT/hist_ensbl/fig_papier/gsubpolmsk.nc").gsubpolmsk

# #MASK GYRE SUBPOLAIRE
# basin = xr.open_mfdataset("Documents/TT/hist_ensbl//python/basins_e*.nc")
# gb = basin['global']

# gsubpolmsk  = np.logical_and(gb, np.abs(gb.y)>=200)
# gsubpolmsk =  np.logical_and(gsubpolmsk, np.abs(gsubpolmsk.x)>=230)
# gsubpolmsk =  np.logical_and(gsubpolmsk, np.abs(gsubpolmsk.y)<=250)
# gsubpolmsk =  np.logical_and(gsubpolmsk, np.abs(gsubpolmsk.x)<=270)

# hist_cor_anom_3000 = xr.open_mfdataset("/home/seguyr/Documents/TT/hist_ensbl/fig_papier/hist_cor/bootstrap/ohc/2D_bootstrap_hist_corr_anom_3000_ohc.nc").__xarray_dataarray_variable__ 
# test = hist_cor_anom_3000.isel(stats=1)

test2 = ((ohc_2D_cor_trends_3000.mean(dim='ensemble') - ohc_2D_cor_trends_1000.mean(dim='ensemble')).__xarray_dataarray_variable__*(10**12))*gsubpolmsk
test2.plot()

#application mask 
#ohc_2D_trends_gsp_1000 = ((ohc_2D_trends_1000*area_oce)*gsubpolmsk).sum(dim=('x', 'y')).__xarray_dataarray_variable__.values
ohc_2D_cor_trends_gsp_1000 = ((ohc_2D_cor_trends_1000*area_oce)*gsubpolmsk).sum(dim=('x', 'y')).__xarray_dataarray_variable__
#ohc_2D_trends_gsp_3000 = ((ohc_2D_trends_3000*area_oce)*gsubpolmsk).sum(dim=('x', 'y')).__xarray_dataarray_variable__.values
ohc_2D_cor_trends_gsp_3000 = ((ohc_2D_cor_trends_3000*area_oce)*gsubpolmsk).sum(dim=('x', 'y')).__xarray_dataarray_variable__
#ohc_2D_trends_gsp_autre = ((ohc_2D_trends_autre*area_oce)*gsubpolmsk).sum(dim=('x', 'y')).__xarray_dataarray_variable__.values
#ohc_2D_cor_trends_gsp_autre = ((ohc_2D_cor_trends_autre*area_oce)*gsubpolmsk).sum(dim=('x', 'y')).__xarray_dataarray_variable__.values


# ohc_2D_cor_trends_gsp_1000 = (ohc_2D_cor_trends_1000*area_oce).sel(x = 252).sel(y = 225).__xarray_dataarray_variable__   
# ohc_2D_cor_trends_gsp_3000 = (ohc_2D_cor_trends_3000*area_oce).sel(x = 252).sel(y = 225).__xarray_dataarray_variable__


boot_1000_ohc = bootstrap(ohc_2D_cor_trends_gsp_1000, n_boot)
boot_1000_ohc = xr.DataArray(boot_1000_ohc,dims=("stats"),coords={"stats": ["lower", "mean", "upper",  "sigma", "p_value"]})

boot_3000_ohc = bootstrap(ohc_2D_cor_trends_gsp_3000, n_boot)
boot_3000_ohc = xr.DataArray(boot_3000_ohc,dims=("stats"),coords={"stats": ["lower", "mean", "upper",  "sigma", "p_value"]})


boot_1000_amoc = bootstrap(AMOC_cor_trends_1000, n_boot)
boot_1000_amoc = xr.DataArray(boot_1000_amoc,dims=("stats"),coords={"stats": ["lower", "mean", "upper",  "sigma", "p_value"]})

boot_3000_amoc = bootstrap(AMOC_cor_trends_3000, n_boot)
boot_3000_amoc = xr.DataArray(boot_3000_amoc,dims=("stats"),coords={"stats": ["lower", "mean", "upper",  "sigma", "p_value"]})


#Moyennes

ohc_cor_mean_1000 = boot_1000_ohc.sel(stats = "mean")
amoc_cor_mean_1000 = boot_1000_amoc.sel(stats = "mean")

ohc_cor_mean_3000 = boot_3000_ohc.sel(stats = "mean")
amoc_cor_mean_3000 = boot_3000_amoc.sel(stats = "mean")


# ohc_cor_mean_1000 = np.mean(ohc_2D_cor_trends_gsp_1000)#/(10**12))
# amoc_cor_mean_1000 = np.mean(AMOC_cor_trends_1000)
# ohc_mean_1000 = np.mean(ohc_2D_trends_gsp_1000)#/(10**12))
# amoc_mean_1000 = np.mean(AMOC_trends_1000)

# ohc_cor_mean_3000 = np.mean(ohc_2D_cor_trends_gsp_3000)#/(10**12))
# amoc_cor_mean_3000 = np.mean(AMOC_cor_trends_3000)
# ohc_mean_3000 = np.mean(ohc_2D_trends_gsp_3000/(10**12))
# amoc_mean_3000 = np.mean(AMOC_trends_3000)

# amoc_mean_mean_1000 = np.mean(AMOC_1900_1950_1000)
# amoc_mean_mean_3000 = np.mean(AMOC_1900_1950_3000)

# amoc_cor_mean_mean_1000 = np.mean(AMOC_cor_1900_1950_1000)
# amoc_cor_mean_mean_3000 = np.mean(AMOC_cor_1900_1950_3000)

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
# x_1000 = np.array([val[()] for val in AMOC_trends_1000])
x_1000 = np.array([val[()] for val in AMOC_cor_trends_1000])

# x_1000 = np.array([val[()] for val in AMOC_1900_1950_1000])
# x_1000 = np.array([val[()] for val in AMOC_cor_1900_1950_1000])
#x_1000 = ohc_2D_trends_gsp_1000/(10**12)

# x_3000 = np.array([val[()] for val in AMOC_trends_3000])
x_3000 = np.array([val[()] for val in AMOC_cor_trends_3000])

# x_3000 = np.array([val[()] for val in AMOC_1900_1950_3000])
# x_3000 = np.array([val[()] for val in AMOC_cor_1900_1950_3000])
#x_3000 = ohc_2D_trends_gsp_3000/(10**12)


# x_autre = np.array([val[()] for val in AMOC_trends_autre])
# x_autre = np.array([val[()] for val in AMOC_cor_trends_autre])
# x_autre = np.array([val[()] for val in AMOC_1900_1950_autre])
# x_autre = np.array([val[()] for val in AMOC_cor_1900_1950_autre])
#x_autre = ohc_2D_trends_gsp_autre/(10**12)


# y_1000 = ohc_2D_trends_gsp_1000/(10**12)
# y_1000 = np.array([val[()] for val in AMOC_cor_trends_1000])
y_1000 = ohc_2D_cor_trends_gsp_1000#/(10**12)


# y_3000 = ohc_2D_trends_gsp_3000/(10**12)
# y_3000 = np.array([val[()] for val in AMOC_cor_trends_3000])
y_3000 = ohc_2D_cor_trends_gsp_3000#/(10**12)


# y_autre = ohc_2D_trends_gsp_autre/(10**12)
# y_autre = np.array([val[()] for val in AMOC_cor_trends_autre])
# y_autre = ohc_2D_cor_trends_gsp_autre/(10**12)


# correlation_coefficient_autre, _ = stats.pearsonr(x_autre, y_autre)
# slope, intercept, r_value, p_value, std_err = stats.linregress(x_autre, y_autre)
# line_autre = slope * x_autre + intercept


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
plt.show()





# plt.figure(figsize=(15, 10))
# plt.scatter(x_3000, y_3000, color='blue', label='hist+3000')
# plt.plot(x_3000, line_3000, color='blue', label='Linear regression hist_dd+3000', linewidth=2)
# plt.plot(amoc_mean_3000, ohc_mean_3000, 'x', color = 'blue', label = 'Ensemble mean', markersize = 30)
# plt.scatter(x_1000, y_1000, color='orange', label='hist+1000')
# plt.plot(x_1000, line_1000, color='orange', label='Linear regression hist_dd+1000', linewidth=2)
# plt.plot(amoc_mean_1000, ohc_mean_1000, 'x', color = 'orange', label = 'Ensemble mean', markersize = 30)
# # plt.scatter(x_autre, y_autre, color='grey', label='Other hist CMIP members')
# # plt.plot(x_autre, line_autre, color='grey', label='Linear regression other hist CMIP members', linewidth=2)
# plt.xlabel('AMOC dedrifted trends (Sv/century)', fontsize=15)
# plt.ylabel('OHC dedrifted trends (ZJ/century)', fontsize=15)
# plt.legend()
# plt.text(min(x_1000)+3, min(y_1000) +5, 
#          f'Correlation coeff: {correlation_coefficient_1000:.2f}', 
#          fontsize=15,color = "orange", fontweight='bold', 
#          verticalalignment='top')
# plt.text(min(x_3000)+3, min(y_3000) +5, 
#          f'Correlation coeff: {correlation_coefficient_3000:.2f}', 
#          fontsize=15,color = "blue",fontweight='bold', 
#          verticalalignment='top')
# # plt.text(min(x_autre)+3, min(y_autre) +5, 
# #          f'Correlation coeff: {correlation_coefficient_autre:.2f}', 
# #          fontsize=15,color = "grey",fontweight='bold', 
# #          verticalalignment='top')
# plt.grid()
# plt.show()


# plt.figure(figsize=(15, 10))
# plt.scatter(x_3000, y_3000, color='blue', label='hist+3000')
# plt.plot(x_3000, line_3000, color='blue', label='Linear regression hist+3000', linewidth=2)
# plt.plot(amoc_mean_mean_3000, ohc_mean_3000, 'x', color = 'blue', label = 'Ensemble mean', markersize = 30)
# plt.scatter(x_1000, y_1000, color='orange', label='hist+1000')
# plt.plot(x_1000, line_1000, color='orange', label='Linear regression hist+1000', linewidth=2)
# plt.plot(amoc_mean_mean_1000, ohc_mean_1000, 'x', color = 'orange', label = 'Ensemble mean', markersize = 30)
# plt.scatter(x_autre, y_autre, color='grey', label='Other hist CMIP members')
# plt.plot(x_autre, line_autre, color='grey', label='Linear regression other hist CMIP members', linewidth=2)
# plt.xlabel('AMOC initial (Sv)', fontsize=15)
# plt.ylabel('OHC trends (ZJ/century)', fontsize=15)
# plt.legend()
# plt.text(min(x_1000)+3, min(y_1000) +5, 
#          f'Correlation coeff: {correlation_coefficient_1000:.2f}', 
#          fontsize=15,color = "orange", fontweight='bold', 
#          verticalalignment='top')
# plt.text(min(x_3000)+3, min(y_3000) +5, 
#          f'Correlation coeff: {correlation_coefficient_3000:.2f}', 
#          fontsize=15,color = "blue",fontweight='bold', 
#          verticalalignment='top')
# plt.text(min(x_autre)+3, min(y_autre) +5, 
#          f'Correlation coeff: {correlation_coefficient_autre:.2f}', 
#          fontsize=15,color = "grey",fontweight='bold', 
#          verticalalignment='top')
# plt.grid()
# plt.show()


# plt.figure(figsize=(15, 10))
# plt.scatter(x_3000, y_3000, color='blue', label='hist+3000')
# plt.plot(x_3000, line_3000, color='blue', label='Linear regression hist_cor+3000', linewidth=2)
# plt.plot(amoc_cor_mean_mean_3000, ohc_cor_mean_3000, 'x', color = 'blue', label = 'Ensemble mean', markersize = 30)
# plt.scatter(x_1000, y_1000, color='orange', label='hist+1000')
# plt.plot(x_1000, line_1000, color='orange', label='Linear regression hist_cor+1000', linewidth=2)
# plt.plot(amoc_cor_mean_mean_1000, ohc_cor_mean_1000, 'x', color = 'orange', label = 'Ensemble mean', markersize = 30)
# # plt.scatter(x_autre, y_autre, color='grey', label='Other hist_cor CMIP members')
# # plt.plot(x_autre, line_autre, color='grey', label='Linear regression other hist_cor CMIP members', linewidth=2)
# plt.xlabel('AMOC cor initial (Sv)', fontsize=15)
# plt.ylabel('OHC cor trends (ZJ/century)', fontsize=15)
# plt.legend()
# plt.text(min(x_1000)+3, min(y_1000) +5, 
#          f'Correlation coeff: {correlation_coefficient_1000:.2f}', 
#          fontsize=15,color = "orange", fontweight='bold', 
#          verticalalignment='top')
# plt.text(min(x_3000)+3, min(y_3000) +20, 
#          f'Correlation coeff: {correlation_coefficient_3000:.2f}', 
#          fontsize=15,color = "blue",fontweight='bold', 
#          verticalalignment='top')
# # plt.text(min(x_autre)+3, min(y_autre) +5, 
# #          f'Correlation coeff: {correlation_coefficient_autre:.2f}', 
# #          fontsize=15,color = "grey",fontweight='bold', 
# #          verticalalignment='top')
# plt.grid()
# plt.show()















