"""
Figure 12: Positinon in CMIP ensemble
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
INTERMEDIATE_DIR = PROJECT_ROOT / "data/data_plot/fig12"



ref_start = 1850
ref_end = 1900
p_start = 1995
p_end = 2014
x_gb, y_gb = 0.7, 0.84
IC = "90%"
var = 'ohc'
unit = 'ZJ'


diff_hist_cor_tot_TM = xr.open_mfdataset(INTERMEDIATE_DIR  / "diff_hist_cor_tot_TM.nc").__xarray_dataarray_variable__
diff_hist_cor_tot_TM_2 = xr.open_mfdataset(INTERMEDIATE_DIR  / "diff_hist_cor_tot_TM_2.nc").__xarray_dataarray_variable__

hist_dd = xr.open_mfdataset(INTERMEDIATE_DIR  / "hist_dd.nc").ohc_full_J_detrended
mean = xr.open_mfdataset(INTERMEDIATE_DIR  / "mean.nc").ohc_full_J_detrended
q05 = xr.open_mfdataset(INTERMEDIATE_DIR  / "q05.nc").ohc_full_J_detrended
q95 = xr.open_mfdataset(INTERMEDIATE_DIR  / "q95.nc").ohc_full_J_detrended

mean_tot = xr.open_mfdataset(INTERMEDIATE_DIR  / "mean_tot.nc").__xarray_dataarray_variable__
q05_tot = xr.open_mfdataset(INTERMEDIATE_DIR  / "q05_tot.nc").__xarray_dataarray_variable__
q95_tot = xr.open_mfdataset(INTERMEDIATE_DIR  / "q95_tot.nc").__xarray_dataarray_variable__

mean_3000 = xr.open_mfdataset(INTERMEDIATE_DIR  / "mean_3000.nc").__xarray_dataarray_variable__
q05_3000 = xr.open_mfdataset(INTERMEDIATE_DIR  / "q05_3000.nc").__xarray_dataarray_variable__
q95_3000 = xr.open_mfdataset(INTERMEDIATE_DIR  / "q95_3000.nc").__xarray_dataarray_variable__
panel_labels = ["a)", "b)"]

time = diff_hist_cor_tot_TM.time + 1850

fig, ax = plt.subplots(
    2, 1,
    figsize=(10, 10),
    sharex=True,
    sharey=False
)

# ============================================================
# SUBPLOT 1
# ============================================================
colors = {"bleu foncé" : '#1b365d', "blue clair":'#a3bce2', "teal foncé":'#005f73', "teal clair":'#9ae3d2', "rouge foncé":'#b2182b', "rouge clair":'#f4a582', "orange foncé":'#e07b00', "orange clair":'#fdd87a' }

ax[0].plot(
    hist_dd.time,
    mean,
    color = "grey",
    label="hist_dd_CMIP6",
)

ax[0].fill_between(
    hist_dd.time,
    q05,
    q95,
    color = "grey",
    alpha=0.2,
)

ax[0].plot(
    hist_dd.time,
    mean_3000,
    color = colors['teal foncé'],
    label="hist_dd+3000",
)

ax[0].fill_between(
    hist_dd.time,
    q05_3000,
    q95_3000,
    color = colors['teal clair']
)


ax[0].plot(
    hist_dd.time,
    mean_tot,
    color = colors['orange foncé'],
    label="hist_dd+1000",
)

ax[0].fill_between(
    hist_dd.time,
    q05_tot,
    q95_tot,
    color = colors['orange clair'],
    alpha = 0.5
)

ax[0].axhline(0, linewidth=0.8, color = "green")


ax[0].set_ylabel("OHC hist_dd (ZJ)", fontsize = 15)
ax[0].legend(loc="upper left", fontsize = 15)
ax[0].tick_params(axis="both", which="major", labelsize=15)

ax[0].text(
        -0.10, 1.02, panel_labels[0],
        transform=ax[0].transAxes,
        ha="right", va="bottom",
        fontsize=20,
        fontweight="bold",
        clip_on=False
    )

# ============================================================
# SUBPLOT 2
# ============================================================

ax[1].plot(
    time,
    diff_hist_cor_tot_TM.sel(stats='mean').values,
    color="red",
    label=f'[Strong drift(13members) - Weak drift(31members)]'
)

ax[1].fill_between(
    time,
    diff_hist_cor_tot_TM.sel(stats='lower').values,
    diff_hist_cor_tot_TM.sel(stats='upper').values,
    color="pink",
    alpha=0.3
)

ax[1].plot(
    time,
    diff_hist_cor_tot_TM_2.sel(stats='mean').values,
    color="blue",
    label=f'[Cold state(28members) - Warm state(16members)]'
)

ax[1].fill_between(
    time,
    diff_hist_cor_tot_TM_2.sel(stats='lower').values,
    diff_hist_cor_tot_TM_2.sel(stats='upper').values,
    color="lightblue",
    alpha=0.3
)

ax[1].axhline(0, linewidth=0.8, color="green")

ax[1].set_xlabel("Simulated years", fontsize = 15)
ax[1].set_ylabel("OHC diff_hist_dd (ZJ)", fontsize = 15)
ax[1].legend(loc="upper left", fontsize = 15)
ax[1].tick_params(axis="both", which="major", labelsize=15)

ax[1].text(
        -0.10, 1.02, panel_labels[1],
        transform=ax[1].transAxes,
        ha="right", va="bottom",
        fontsize=20,
        fontweight="bold",
        clip_on=False
    )
# ============================================================
# Mise en page
# ============================================================

plt.tight_layout()


plt.savefig(FIG_DIR / "figure_12.pdf", bbox_inches="tight")
plt.savefig(FIG_DIR / "figure_12.png", bbox_inches="tight")

