"""
Figure 11: Influence of ensemble size 
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
INTERMEDIATE_DIR = PROJECT_ROOT / "data/data_plot/fig11"

m_OHC_dd_boot_both = xr.open_mfdataset(INTERMEDIATE_DIR / "m_OHC_dd_boot_both.nc").__xarray_dataarray_variable__


#Tracé

# Axe croissant : 10 à 29 membres
x = np.arange(4, 30, 2)

# Les données sont stockées dans l'ordre 29 -> 10
# donc on inverse leur ordre
mean = m_OHC_dd_boot_both.sel(stats="mean").values[::-1]
lower = m_OHC_dd_boot_both.sel(stats="lower").values[::-1]
upper = m_OHC_dd_boot_both.sel(stats="upper").values[::-1]

# IC à 90 %
ic_lower = lower
ic_upper = upper

fig, ax = plt.subplots(figsize=(8, 5))

# Moyenne
ax.plot(x, mean, "-", color = 'red', linewidth=2, label="Mean")

# IC 90 %
ax.fill_between(
    x,
    ic_lower,
    ic_upper,
    color = "red",
    alpha=0.3,
    label="90% IC"
)

ax.set_xlabel("Number of members")
ax.set_ylabel("OHC [hist_dd+3000 - hist_dd+1000] (in ZJ)")
ax.set_xticks(x)

ax.axhline(0, color='green', linewidth=3)

ax.legend(loc = "lower right")
ax.grid(alpha=0.3)

plt.tight_layout()


plt.savefig(FIG_DIR / "figure_11.pdf", bbox_inches="tight")
plt.savefig(FIG_DIR / "figure_11.png", bbox_inches="tight")
