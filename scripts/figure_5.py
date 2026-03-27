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
# Colors
# -----------------------------------------------------------------------------
colors = {
    "orange foncé": "#d95f02",
    "orange clair": "#fdd0a2",
    "teal foncé": "#1b9e77",
    "teal clair": "#a6dba0",
    "rouge foncé": "#b2182b",
    "rouge clair": "#f4a6b3",
    "bleu foncé": "#2166ac",
    "blue clair": "#92c5de",
}


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
def integrate_global_ohc(arr, area_oce, scale=1.0):
    """
    Convert a 2D OHC field (J m^-2) into a global OHC time series (scaled).

    Parameters
    ----------
    arr : xr.DataArray
        Input array with dimensions including x and y.
    area_oce : xr.DataArray
        Ocean cell area (m²).
    scale : float
        Scaling factor applied after integration.

    Returns
    -------
    xr.DataArray
        Integrated global OHC time series.
    """
    spatial_dims = [d for d in arr.dims if d not in {"ensemble", "time", "stats"}]
    return (arr * area_oce).sum(dim=spatial_dims) / scale


def plot_ci(ax, time, arr, line_color, fill_color, alpha=0.8, lw=2.5):
    """
    Plot bootstrap mean and confidence interval from a DataArray with
    a 'stats' dimension.
    """
    lower, mean, upper = get_stats(arr)

    ax.plot(time, mean.values, color=line_color, linewidth=lw)
    ax.fill_between(
        time,
        lower.values,
        upper.values,
        color=fill_color,
        alpha=alpha,
    )


def main():
    """Generate Figure 5."""

    # -------------------------------------------------------------------------
    # Load data
    # -------------------------------------------------------------------------
    area_oce = load_area_ocean()

    (
        OHC_2D_hist_tot,
        OHC_2D_hist_3000,
        OHC_2D_pic_tot,
        OHC_2D_pic_3000,
    ) = load_ohc_2d_ensembles()

    # -------------------------------------------------------------------------
    # Global OHC time series for raw ensembles
    # -------------------------------------------------------------------------
    pic_1000 = integrate_global_ohc(OHC_2D_pic_tot, area_oce, scale=ECHELLE_OHC)
    pic_3000 = integrate_global_ohc(OHC_2D_pic_3000, area_oce, scale=ECHELLE_OHC)

    hist_1000 = integrate_global_ohc(OHC_2D_hist_tot, area_oce, scale=ECHELLE_OHC)
    hist_3000 = integrate_global_ohc(OHC_2D_hist_3000, area_oce, scale=ECHELLE_OHC)

    # Anomalies relative to first 50 years (1850-1899)
    pic_1000 = anomalies(pic_1000)
    pic_3000 = anomalies(pic_3000)
    hist_1000 = anomalies(hist_1000)
    hist_3000 = anomalies(hist_3000)

    # Bootstrap raw time series
    pic_1000 = boot(pic_1000, n_boot=N_BOOT)
    pic_3000 = boot(pic_3000, n_boot=N_BOOT)
    hist_1000 = boot(hist_1000, n_boot=N_BOOT)
    hist_3000 = boot(hist_3000, n_boot=N_BOOT)

    # -------------------------------------------------------------------------
    # Dedrifted historical response: hist - matched pic
    # Keep 2D fields first, then integrate globally
    # -------------------------------------------------------------------------
    OHC_dd_1000 = time_matching(OHC_2D_hist_tot, OHC_2D_pic_tot)
    OHC_dd_3000 = time_matching(OHC_2D_hist_3000, OHC_2D_pic_3000)

    OHC_dd_1000 = anomalies(OHC_dd_1000)
    OHC_dd_3000 = anomalies(OHC_dd_3000)

    ohc_dd_1000 = integrate_global_ohc(OHC_dd_1000, area_oce, scale=ECHELLE_OHC)
    ohc_dd_3000 = integrate_global_ohc(OHC_dd_3000, area_oce, scale=ECHELLE_OHC)

    # Mean gain over final 20 years
    m_ohc_dd_1000 = gain(ohc_dd_1000)
    m_ohc_dd_3000 = gain(ohc_dd_3000)

    # Bootstrap dedrifted time series and gain differences
    diff_hist_cor = boot_diff(ohc_dd_1000, ohc_dd_3000, n_boot=N_BOOT)
    m_diff_hist_cor = boot_diff(m_ohc_dd_1000, m_ohc_dd_3000, n_boot=N_BOOT)

    # Extract scalar stats for the late-period mean difference
    rm_diff_low, rm_diff_mean, rm_diff_up = get_scalar_stats(m_diff_hist_cor)

    # -------------------------------------------------------------------------
    # Figure
    # -------------------------------------------------------------------------
    fig, axes = plt.subplots(2, 1, figsize=(15, 15), sharex=True)
    ax1, ax2 = axes

    panel_labels = ["a)", "b)"]

    # ============================================================
    # SUBPLOT 1: raw historical and piControl-derived responses
    # ============================================================
    time = pic_1000.time.values + REF_START

    plot_ci(
        ax1, time, pic_1000,
        line_color=colors["rouge foncé"],
        fill_color=colors["rouge clair"],
    )
    plot_ci(
        ax1, time, pic_3000,
        line_color=colors["bleu foncé"],
        fill_color=colors["blue clair"],
    )
    plot_ci(
        ax1, time, hist_1000,
        line_color=colors["orange foncé"],
        fill_color=colors["orange clair"],
    )
    plot_ci(
        ax1, time, hist_3000,
        line_color=colors["teal foncé"],
        fill_color=colors["teal clair"],
    )

    ax1.set_ylabel(f"OHC anomaly relative to {REF_START}-{REF_END} ({UNIT})", fontsize=20)
    ax1.set_title("CNRM-CM6.1 working ensembles", fontsize=25, fontweight="bold")
    ax1.tick_params(axis="both", labelsize=18)
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
    # SUBPLOT 2: difference of dedrifted historical responses
    # ============================================================
    diff_low, diff_mean, diff_up = get_stats(diff_hist_cor)
    time = diff_mean.time.values + REF_START

    ax2.plot(time, diff_mean.values, color="grey", linewidth=2.5)
    ax2.fill_between(
        time,
        diff_low.values,
        diff_up.values,
        color="red",
        alpha=0.3,
    )

    # Late-period mean difference
    ax2.hlines(
        rm_diff_mean,
        PERIOD_GAIN_START,
        PERIOD_GAIN_END,
        colors="purple",
        linewidth=3,
    )
    ax2.errorbar(
        x=0.5 * (PERIOD_GAIN_START + PERIOD_GAIN_END),
        y=rm_diff_mean,
        yerr=[[rm_diff_mean - rm_diff_low], [rm_diff_up - rm_diff_mean]],
        fmt="o",
        color="purple",
        capsize=5,
    )

    ax2.set_ylabel(f"OHC anomaly relative to {REF_START}-{REF_END} ({UNIT})", fontsize=20)
    ax2.set_xlabel("Year", fontsize=20)
    ax2.tick_params(axis="both", labelsize=18)
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

    # ============================================================
    # Legend
    # ============================================================
    legend_handles = [
        Patch(facecolor=colors["orange foncé"], edgecolor="none", label="Hist_dd+1000"),
        Patch(facecolor=colors["teal foncé"], edgecolor="none", label="Hist_dd+3000"),
        Patch(facecolor=colors["rouge foncé"], edgecolor="none", label="Pic+1000"),
        Patch(facecolor=colors["bleu foncé"], edgecolor="none", label="Pic+3000"),
        Patch(facecolor=colors["rouge clair"], edgecolor="none", label="Hist_dd+3000 - Hist_dd+1000"),
    ]

    fig.legend(
        handles=legend_handles,
        loc="upper center",
        ncol=5,
        frameon=False,
        fontsize=18,
        bbox_to_anchor=(0.5, -0.01),
    )

    plt.tight_layout(h_pad=3)

    # Save
    plt.savefig(FIG_DIR / "figure5.pdf", bbox_inches="tight")
    plt.savefig(FIG_DIR / "figure5.png", dpi=300, bbox_inches="tight")

    plt.show()


if __name__ == "__main__":
    main()
