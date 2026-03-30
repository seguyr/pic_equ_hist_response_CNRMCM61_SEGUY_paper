"""
Figure 2: piControl equilibrium diagnostics
This script plots:
- TOA flux over the last 1500 simulated years,
- rolling temperature trends (SAT, SST, DOT),
- rolling MOC trends (AMOC, PMOC, GMOC, SOMOC).
"""

from pathlib import Path
import sys
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

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

# -----------------------------------------------------------------------------
# Imports from project utilities
# -----------------------------------------------------------------------------
from functions.utils import (  # noqa: E402
    TAS_PIC_FILE,
    TOS_PIC_FILE,
    MOC_PIC_FILE,
    DOT_PIC_FILE,
    TOA_PIC_FILE,
    rolling_trend_np,
    rolling_mean_np,
)

# -----------------------------------------------------------------------------
# Matplotlib style
# -----------------------------------------------------------------------------
plt.rc("font", size=12.5)
plt.rcParams["figure.dpi"] = 300

# -----------------------------------------------------------------------------
# Load data
# -----------------------------------------------------------------------------
SAT = xr.open_dataset(TAS_PIC_FILE)
vdata_sat = SAT[list(SAT.data_vars)[0]]

SST = xr.open_dataset(TOS_PIC_FILE)
vdata_sst = SST[list(SST.data_vars)[0]]

MOC = xr.open_dataset(MOC_PIC_FILE)
vdata_moc = MOC

DOT = xr.open_dataset(DOT_PIC_FILE)
vdata_dot = DOT[list(DOT.data_vars)[0]]

TOA = xr.open_dataset(TOA_PIC_FILE)
vdata_toa = TOA[list(TOA.data_vars)[0]]
vdata_toa = vdata_toa.assign_coords(year=np.arange(vdata_toa.sizes["year"]))

# -----------------------------------------------------------------------------
# Rolling diagnostics
# -----------------------------------------------------------------------------
slopes, years = rolling_trend_np(vdata_sat.values, window=1500)
slopes_100y_sat = slopes * 100

slopes, years = rolling_trend_np(vdata_sst.values, window=1500)
slopes_100y_sst = slopes * 100

mean_1500y_toa = rolling_mean_np(vdata_toa.values, window=1500)
slopes, years = rolling_trend_np(vdata_toa.values, window=1500)
slopes_100y_toa = slopes * 100

slopes, years = rolling_trend_np(vdata_dot.values, window=1500)
slopes_100y_dot = slopes * 100

slopes, years = rolling_trend_np(vdata_moc["AMOC"].values, window=1500)
slopes_1000y_amoc = slopes * 1000

slopes, years = rolling_trend_np(vdata_moc["PMOC"].values, window=1500)
slopes_1000y_pmoc = slopes * 1000

slopes, years = rolling_trend_np(vdata_moc["SOMOC"].values, window=1500)
slopes_1000y_somoc = slopes * 1000

slopes, years = rolling_trend_np(vdata_moc["GMOC"].values, window=1500)
slopes_1000y_gmoc = slopes * 1000

# -----------------------------------------------------------------------------
# Figure
# -----------------------------------------------------------------------------
fig, axs = plt.subplots(3, 1, figsize=(15, 15))
panel_labels = ["a)", "b)", "c)"]

# -------------------------------------------------------------------------
# Panel a: TOA
# -------------------------------------------------------------------------
axs[0].plot(
    np.arange(1500, 3000),
    vdata_toa.isel(year=slice(1500, 3000)),
    color="grey",
)
axs[0].plot(
    np.arange(1499, 3001),
    mean_1500y_toa,
    color="black",
    label="TOA rolling mean",
)
axs[0].hlines(0.1, 1500, 3000, color="red", linewidth=5)
axs[0].hlines(0.1, 1500, 3000, color="green", linewidth=2, label="< 0.1 W m$^{-2}$")
axs[0].set_ylabel("TOA flux (W m$^{-2}$)", fontsize=20, fontweight="bold")
axs[0].tick_params(axis="both", which="major", labelsize=20)
axs[0].set_title("Radiative balance", fontsize=25, fontweight="bold")
axs[0].grid(True)
axs[0].legend(loc="lower right", fontsize=15)

axs[0].text(
    -0.08, 1.02, panel_labels[0],
    transform=axs[0].transAxes,
    ha="right", va="bottom",
    fontsize=22, fontweight="bold",
    clip_on=False,
)

# -------------------------------------------------------------------------
# Panel b: temperatures
# -------------------------------------------------------------------------
axs[1].plot(years, slopes_100y_sat, color="black", label="SAT")
axs[1].plot(years, slopes_100y_sst, color="grey", label="SST")
axs[1].plot(years, slopes_100y_dot, "--", color="grey", label="DOT")
axs[1].hlines(0.05, 1500, 3000, color="red", linewidth=5)
axs[1].hlines(0.1, 1500, 3000, color="green", linewidth=2, label="< 0.1 °C century$^{-1}$")
axs[1].set_ylabel("1500y trends (°C century$^{-1}$)", fontsize=20, fontweight="bold")
axs[1].tick_params(axis="both", which="major", labelsize=20)
axs[1].set_title("Temperature", fontsize=25, fontweight="bold")
axs[1].grid(True)
axs[1].legend(fontsize=15)

axs[1].text(
    -0.08, 1.02, panel_labels[1],
    transform=axs[1].transAxes,
    ha="right", va="bottom",
    fontsize=22, fontweight="bold",
    clip_on=False,
)

# -------------------------------------------------------------------------
# Panel c: MOC
# -------------------------------------------------------------------------
axs[2].plot(years, slopes_1000y_amoc, color="black", label="AMOC")
axs[2].plot(years, slopes_1000y_pmoc, "--", color="black", label="PMOC")
axs[2].plot(years, slopes_1000y_gmoc, color="black", linewidth=3, label="GMOC")
axs[2].plot(years, slopes_1000y_somoc, color="grey", label="SOMOC")
axs[2].hlines(1, 1500, 3000, color="green", linewidth=2)
axs[2].hlines(-1, 1500, 3000, color="green", linewidth=2)
axs[2].fill_between(years, 1, -1, color="lightgreen", alpha=0.5, label="< 1 Sv millennial$^{-1}$")
axs[2].set_xlabel("Simulated years", fontsize=20)
axs[2].tick_params(axis="both", which="major", labelsize=20)
axs[2].set_ylabel("1500y trends (Sv millennial$^{-1}$)", fontsize=20, fontweight="bold")
axs[2].set_title("MOC", fontsize=25, fontweight="bold")
axs[2].grid(True)
axs[2].legend(fontsize=15)

axs[2].text(
    -0.08, 1.02, panel_labels[2],
    transform=axs[2].transAxes,
    ha="right", va="bottom",
    fontsize=22, fontweight="bold",
    clip_on=False,
)

plt.tight_layout()

plt.savefig(FIG_DIR / "figure_2.pdf", bbox_inches="tight")














