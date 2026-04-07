
"""
Figure 8: obs VS model comparison for OHC 
====================================================

"""

from pathlib import Path
import sys
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch

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
)

ds = xr.open_mfdataset(INTERMEDIATE_DIR / "OHC_prepared_3layers_M2.nc")
ds_bar = xr.open_dataset(INTERMEDIATE_DIR / "barplot_values_yerr.nc")
yerr   = ds_bar["yerr"].values
values = ds_bar["values"].values

P_START = 1989
P_END = 2014
REF_START = 1958
REF_END = 1983

t_m = np.linspace(1850,2014,165)
t_o = ds["time_obs"].values
layers = ds["layer"].values

layer_labels = ["0–300 m", "300–2000 m", "2000–bottom"]
bar_labels = ["Hist_dd+1000", "Hist_dd+3000", "OHC Reference products"]

n_layers = len(layer_labels)
n_bars   = len(bar_labels)
xpos = np.arange(n_bars)


bar_colors = [
    COLORS["orange_dark"], 
    COLORS["teal_dark"],   
    COLORS["red_dark"],   
]

panel_labels = ["a)", "b)", "c)"]


# ============================================================
# FIGURE 1x2
# ============================================================


fig = plt.figure(figsize=(18, 12), constrained_layout=True)
gs = gridspec.GridSpec(n_layers, 2, figure=fig, width_ratios=[3.8, 1.6])

axes_ts  = [fig.add_subplot(gs[i, 0]) for i in range(n_layers)]
axes_bar = [fig.add_subplot(gs[i, 1]) for i in range(n_layers)]

# ============================================================
# TITRES DE COLONNES
# ============================================================
axes_ts[0].set_title("OHC (ZJ)", fontsize=20, fontweight="bold", pad=18)
axes_bar[0].set_title("OHC response over 1958–2014 (ZJ)", fontsize=20, fontweight="bold", pad=18)

cap_size = 10
err_lw   = 2.2
cap_thk  = 2.2

for i in range(n_layers):
    # =========================
    # (A) TIME SERIES
    # =========================
    ax = axes_ts[i]
    layer = layers[i]

    b1000 = ds["model_dd_1000_boot"].sel(layer=layer)
    b3000 = ds["model_dd_3000_boot"].sel(layer=layer)

    ax.plot(
        t_m, b1000.sel(stats="mean").values,
        color=COLORS["orange_dark"], lw=2, label="Hist_dd+1000"
    )
    ax.fill_between(
        t_m,
        b1000.sel(stats="lower").values,
        b1000.sel(stats="upper").values,
        color=COLORS["orange_light"], alpha=0.5
    )

    ax.plot(
        t_m, b3000.sel(stats="mean").values,
        color=COLORS["teal_dark"], lw=2, label="Hist_dd+3000"
    )
    ax.fill_between(
        t_m,
        b3000.sel(stats="lower").values,
        b3000.sel(stats="upper").values,
        color=COLORS["teal_light"], alpha=0.5
    )

    om = ds["obs_mean"].sel(layer=layer)
    oi = ds["obs_ic"].sel(layer=layer)
    if np.isfinite(om.values).any():
        ax.plot(t_o, om.values, color=COLORS["red_dark"], lw=2, label="OHC Reference products")
        ax.fill_between(
            t_o, (om - oi).values, (om + oi).values,
            color=COLORS["red_light"], alpha=0.6
        )

    ax.grid(True, alpha=0.25)
    ax.tick_params(axis="both", which="major", labelsize=16)

    if "p_start" in ds.attrs and "p_end" in ds.attrs:
        ax.axvspan(ds.attrs["p_start"], ds.attrs["p_end"],
                   color="0.85", alpha=0.5, zorder=0)

    if i == n_layers - 1:
        ax.set_xlabel("Year", fontsize=20)
    else:
        ax.set_xlabel("")

    ax.text(
        -0.10, 1.02, panel_labels[i],
        transform=ax.transAxes,
        ha="right", va="bottom",
        fontsize=20,
        fontweight="bold",
        clip_on=False
    )
    #ax.axvspan(REF_START, REF_END, color='grey', alpha=0.2)
    #ax.axvspan(P_START, P_END, color='grey', alpha=0.2)
    
    # =========================
    # (B) BARRES + ERROR BARS
    # =========================
    axb = axes_bar[i]

    for sp in ["top", "left", "bottom"]:
        axb.spines[sp].set_visible(False)

    axb.spines["right"].set_visible(True)
    axb.spines["right"].set_linewidth(1.5)

    axb.set_xticks([])
    axb.set_xlabel("")

    axb.yaxis.set_ticks_position("right")
    axb.yaxis.set_label_position("right")
    axb.tick_params(axis="y", labelsize=16)

    axb.grid(True, axis="y", alpha=0.25)

    bar_w = 0.52
    ymins, ymaxs = [], []

    for j in range(n_bars):
        val = values[i, j]
        if np.isnan(val):
            continue

        # --- BARRE
        axb.bar(
            xpos[j], val,
            width=bar_w,
            color=bar_colors[j],
            edgecolor="none",
            alpha=0.60,
            zorder=2
        )

        # --- ERROR BARS asymétriques
        err = yerr[:, i, j]
        has_ci = np.all(np.isfinite(err))

        if has_ci:
            axb.errorbar(
                xpos[j], val,
                yerr=[[err[0]], [err[1]]],
                fmt="none",
                ecolor=bar_colors[j],
                elinewidth=err_lw,
                capsize=cap_size,
                capthick=cap_thk,
                zorder=3
            )
            ymins.append(min(val - err[0], 0))
            ymaxs.append(max(val + err[1], 0))
        else:
            ymins.append(min(val, 0))
            ymaxs.append(max(val, 0))

    if ymins:
        ymin, ymax = min(ymins), max(ymaxs)
        span = ymax - ymin if ymax > ymin else max(abs(ymax), 1.0)
        axb.set_ylim(ymin - 0.10 * span, ymax + 0.15 * span)

    axb.set_xlim(-0.7, n_bars - 1 + 0.7)
    axb.axhline(0, color="black", lw=1, alpha=0.5, zorder=1)

# ============================================================
# TITRES VERTICAUX PAR LIGNE
# ============================================================
for i, label in enumerate(layer_labels):
    axes_ts[i].text(
        -0.14, 0.5, label,
        rotation=90,
        va="center", ha="center",
        transform=axes_ts[i].transAxes,
        fontsize=20,
        fontweight="bold"
    )

# ============================================================
# LÉGENDE UNIQUE (barres)
# ============================================================
bar_handles = [
    Patch(facecolor=bar_colors[k], edgecolor="none", label=bar_labels[k])
    for k in range(n_bars)
]

fig.legend(
    handles=bar_handles,
    loc="upper center",
    ncol=3,
    frameon=False,
    fontsize=25,
    bbox_to_anchor=(0.5, -0.01)
)

plt.savefig(FIG_DIR / "figure_8.pdf", bbox_inches="tight")
plt.savefig(FIG_DIR / "figure_8.png", bbox_inches="tight")