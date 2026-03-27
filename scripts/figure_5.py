"""
Figure 5: Historical time series response comparison
====================================================

This script plots:
- pic_tot time series in anomalies relative to 1850-1900,
- pic_3000 time series in anomalies relative to 1850-1900,
- hist_tot time series in anomalies relative to 1850-1900,
- hist_3000 time series in anomalies relative to 1850-1900,
- the difference between historical OHC responses after dedrifting by time matching.
"""

from pathlib import Path
import sys
import numpy as np
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
from functions.utils import (
    load_ohc_2d_ensembles,
    time_matching,
    anomalies,
    gain,
    boot,
    boot_diff,
    get_stats,
    get_scalar_stats,
    COLORS
)

# -----------------------------------------------------------------------------
# Physical constants
# -----------------------------------------------------------------------------
ECHELLE_OHC = 1e21  # convert J to ZJ
UNIT = "ZJ"

REF_START = 1850
REF_END = 1899

PERIOD_GAIN_START = 1995
PERIOD_GAIN_END = 2014

TIME_HIST = np.arange(165) + REF_START

# -----------------------------------------------------------------------------
# Calcul
# -----------------------------------------------------------------------------

hist_1000 = load_integrated_ohc("hist_tot") / ECHELLE_OHC
hist_3000 = load_integrated_ohc("hist_3000") / ECHELLE_OHC
pic_1000 = load_integrated_ohc("pic_tot") / ECHELLE_OHC
pic_3000 = load_integrated_ohc("pic_3000") / ECHELLE_OHC

# Dedrift by time matching
OHC_dd_1000 = time_matching(hist_1000, pic_1000)
OHC_dd_3000 = time_matching(hist_3000, pic_3000)

# Anomalies relative to first 50 years
OHC_dd_1000 = anomalies(OHC_dd_1000)
OHC_dd_3000 = anomalies(OHC_dd_3000)

# Mean gain over final 20 years
m_OHC_dd_1000 = gain(OHC_dd_1000)
m_OHC_dd_3000 = gain(OHC_dd_3000)

# Bootstrap
hist_cor_anom_1000 = boot(OHC_dd_1000)
hist_cor_anom_3000 = boot(OHC_dd_3000)
diff_hist_cor = boot_diff(OHC_dd_1000, OHC_dd_3000)

m_hist_cor_1000 = boot(m_OHC_dd_1000)
m_hist_cor_3000 = boot(m_OHC_dd_3000)
m_diff_hist_cor = boot_diff(m_OHC_dd_1000, m_OHC_dd_3000)

pic_1000_low, pic_1000_mean, pic_1000_up = get_stats(pic_1000)
pic_3000_low, pic_3000_mean, pic_3000_up = get_stats(pic_3000)

hist_1000_low, hist_1000_mean, hist_1000_up = get_stats(hist_1000)
hist_3000_low, hist_3000_mean, hist_3000_up = get_stats(hist_3000)

diff_low, diff_mean, diff_up = get_stats(diff_hist_cor)
rm_diff_low, rm_diff_mean, rm_diff_up = get_scalar_stats(m_diff_hist_cor)



# -----------------------------------------------------------------------------
# Tracé
# -----------------------------------------------------------------------------

fig, axes = plt.subplots(2, 1, figsize=(15, 15), sharex=True)
ax1, ax2 = axes

bar_colors = [
    COLORS["orange foncé"],
    COLORS["teal foncé"],
    COLORS["rouge foncé"],
    COLORS["bleu foncé"],
    COLORS["rouge clair"],
]
bar_labels = [
    "Hist_dd+1000",
    "Hist_dd+3000",
    "Pic+1000",
    "PiC+3000",
    "Hist_dd+3000 - Hist_dd+1000",
]
panel_labels = ["a)", "b)"]

# ============================================================
# SUBPLOT 1
# ============================================================
time = pic_1000_mean.time + REF_START

ax1.plot(time, pic_1000_mean.values, color=COLORS["rouge foncé"])
ax1.fill_between(
    time,
    pic_1000_low.values,
    pic_1000_up.values,
    color=COLORS["rouge clair"],
    alpha=0.8,
)

ax1.plot(time, pic_3000_mean.values, color=COLORS["bleu foncé"])
ax1.fill_between(
    time,
    pic_3000_low.values,
    pic_3000_up.values,
    color=COLORS["blue clair"],
    alpha=0.8,
)


ax1.plot(time, hist_1000_mean.values, color=COLORS["orange foncé"])
ax1.fill_between(
    time,
    hist_1000_low.values,
    hist_1000_up.values,
    color=COLORS["orange clair"],
    alpha=0.8,
)

ax1.plot(time, hist_3000_mean.values, color=COLORS["teal foncé"])
ax1.fill_between(
    time,
    hist_3000_low.values,
    hist_3000_up.values,
    color=COLORS["teal clair"],
    alpha=0.8,
)

ax1.set_ylabel(f"OHC / {REF_START}-{REF_END} ({UNIT})", fontsize=20)
ax1.set_title("CNRM-CM6.1 working ensembles", fontsize=25, fontweight="bold")
ax1.tick_params(axis="both", labelsize=20)
ax1.axvspan(REF_START, REF_END, color="grey", alpha=0.2)
ax1.grid()

ax1.text(
    -0.08, 1.02, panel_labels[0],
    transform=ax1.transAxes,
    ha="right", va="bottom",
    fontsize=22,
    fontweight="bold",
    clip_on=False,
)

# ============================================================
# SUBPLOT 2
# ============================================================
time = diff_mean.time + REF_START

ax2.plot(time, diff_mean.values, color="grey")
ax2.fill_between(
    time,
    diff_low.values,
    diff_up.values,
    color="red",
    alpha=0.3,
)

ax2.hlines(rm_diff_mean, 1995, 2014, colors="purple", linewidth=3)
ax2.errorbar(
    x=2005,
    y=rm_diff_mean,
    yerr=[[rm_diff_mean - rm_diff_low], [rm_diff_up - rm_diff_mean]],
    fmt="o",
    color="purple",
    capsize=5,
)

ax2.set_ylabel(f"OHC / {REF_START}-{REF_END} ({UNIT})", fontsize=20)
ax2.set_xlabel("Year", fontsize=20)
ax2.tick_params(axis="both", labelsize=20)
ax2.axvspan(REF_START, REF_END, color="grey", alpha=0.2)
ax2.axvspan(PERIOD_GAIN_START, PERIOD_GAIN_END, color="grey", alpha=0.2)
ax2.axhline(0, color="green", linewidth=3)
ax2.grid()
ax2.set_title("Hist_dd+3000 - Hist_dd+1000", fontsize=25, fontweight="bold")

ax2.text(
    -0.08, 1.02, panel_labels[1],
    transform=ax2.transAxes,
    ha="right", va="bottom",
    fontsize=22,
    fontweight="bold",
    clip_on=False,
)

bar_handles = [
    Patch(facecolor=bar_colors[k], edgecolor="none", label=bar_labels[k])
    for k in range(len(bar_labels))
]

fig.legend(
    handles=bar_handles,
    loc="upper center",
    ncol=5,
    frameon=False,
    fontsize=20,
    bbox_to_anchor=(0.5, -0.01),
)

plt.tight_layout(h_pad=3)
plt.show()
