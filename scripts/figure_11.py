"""
Figure 11: boxplot ditrib AMOC VS OHC
====================================================

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


boot_amoc_1000 = xr.open_mfdataset(INTERMEDIATE_DIR / "boot_amoc_1000_boxplot.nc").__xarray_dataarray_variable__
boot_amoc_3000 = xr.open_mfdataset(INTERMEDIATE_DIR / "boot_amoc_3000_boxplot.nc").__xarray_dataarray_variable__
boot_ohc_1000 = xr.open_mfdataset(INTERMEDIATE_DIR / "boot_ohc_1000_boxplot.nc").__xarray_dataarray_variable__
boot_ohc_3000 = xr.open_mfdataset(INTERMEDIATE_DIR / "boot_ohc_3000_boxplot.nc").__xarray_dataarray_variable__

ts_1000 = xr.open_mfdataset(INTERMEDIATE_DIR / "ts_1000_amoc_boxplot.nc").__xarray_dataarray_variable__
ts_3000 = xr.open_mfdataset(INTERMEDIATE_DIR / "ts_3000_amoc_boxplot.nc").__xarray_dataarray_variable__
ts_1000_ohc = xr.open_mfdataset(INTERMEDIATE_DIR / "ts_1000_ohc_boxplot.nc").__xarray_dataarray_variable__
ts_3000_ohc = xr.open_mfdataset(INTERMEDIATE_DIR / "ts_3000_ohc_boxplot.nc").__xarray_dataarray_variable__


mean_1000 = boot_amoc_1000.sel(stats="mean").mean()
low_1000  = boot_amoc_1000.sel(stats="lower").mean()
up_1000   = boot_amoc_1000.sel(stats="upper").mean()

mean_3000 = boot_amoc_3000.sel(stats="mean").mean()
low_3000  = boot_amoc_3000.sel(stats="lower").mean()
up_3000   = boot_amoc_3000.sel(stats="upper").mean()

mean_1000_ohc = boot_ohc_1000.sel(stats="mean").mean()
low_1000_ohc  = boot_ohc_1000.sel(stats="lower").mean()
up_1000_ohc   = boot_ohc_1000.sel(stats="upper").mean()

mean_3000_ohc = boot_ohc_3000.sel(stats="mean").mean()
low_3000_ohc  = boot_ohc_3000.sel(stats="lower").mean()
up_3000_ohc   = boot_ohc_3000.sel(stats="upper").mean()



fig, axs = plt.subplots(1, 2, figsize=(10,6), sharey=False)

# ======================
# SUBPLOT 2 : AMOC
# ======================
ax = axs[1]

ax.boxplot(
    [ts_1000, ts_3000],
    positions=[1,2],
    widths=0.5,
    whis=[5,95],
    showfliers=False,
    boxprops=dict(color="grey"),
    medianprops=dict(color="grey"),
    whiskerprops=dict(color="grey"),
    capprops=dict(color="grey"),
)

ax.errorbar(
    1,
    mean_1000,
    yerr=[[mean_1000-low_1000],[up_1000-mean_1000]],
    fmt='o',
    color='orange',
    capsize=5,
    label="piC+1000 bootstrap IC90"
)

ax.errorbar(
    2,
    mean_3000,
    yerr=[[mean_3000-low_3000],[up_3000-mean_3000]],
    fmt='o',
    color='blue',
    capsize=5,
    label="piC+3000 bootstrap IC90"
)

ax.set_xticks([1,2])
ax.set_xticklabels(["piC+1000","piC+3000"], fontsize = 15)
ax.set_ylabel("AMOC (Sv)", fontsize = 15)
ax.tick_params(axis='both', labelsize=15)
ax.grid(alpha=0.3)

# ======================
# SUBPLOT 1 : OHC
# ======================
ax = axs[0]

ax.boxplot(
    [ts_1000_ohc, ts_3000_ohc],
    positions=[1,2],
    widths=0.5,
    whis=[5,95],
    showfliers=False,
    boxprops=dict(color="black"),
    medianprops=dict(color="black"),
    whiskerprops=dict(color="black"),
    capprops=dict(color="black"),
)

ax.errorbar(
    1,
    mean_1000_ohc,
    yerr=[[mean_1000_ohc-low_1000_ohc],[up_1000_ohc-mean_1000_ohc]],
    fmt='o',
    color='orange',
    capsize=5,
)

ax.errorbar(
    2,
    mean_3000_ohc,
    yerr=[[mean_3000_ohc-low_3000_ohc],[up_3000_ohc-mean_3000_ohc]],
    fmt='o',
    color='blue',
    capsize=5,
)

ax.set_xticks([1,2])
ax.set_xticklabels(["piC+1000","piC+3000"], fontsize = 15)
ax.set_ylabel("OHC (ZJ)", fontsize = 15)
ax.tick_params(axis='both', labelsize=15)

ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig(FIG_DIR / "figure_11.pdf", bbox_inches="tight")
plt.savefig(FIG_DIR / "figure_11.png", bbox_inches="tight")