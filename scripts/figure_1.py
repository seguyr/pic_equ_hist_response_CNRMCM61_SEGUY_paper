"""
Figure 1: Experiment protocol
=============================

This script plots:
- the reference piControl global ocean heat content (OHC),
- the reference AMOC time series,
- the historical ensemble members branched after 1000 years spin-up,
- the historical ensemble members branched after 3000 years spin-up.

"""

from pathlib import Path
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from matplotlib.lines import Line2D
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

# -----------------------------------------------------------------------------
# Imports from project utilities
# -----------------------------------------------------------------------------
from functions.utils import (
    load_pic_ohc_1D,
    load_pic_amoc,
    load_integrated_ohc,
    load_branching_years,
)

# -----------------------------------------------------------------------------
# Physical constants
# -----------------------------------------------------------------------------
ECHELLE_OHC = 1e21  # ZJ

# -----------------------------------------------------------------------------
# Time axes
# -----------------------------------------------------------------------------
TIME_HIST = np.arange(165)   # 0 to 164

# -----------------------------------------------------------------------------
# Shaded branching windows
# -----------------------------------------------------------------------------
YINI_1000 = 1883 - 1850
YEND_1000 = 2831 - 1850
YINI_3000 = 4462 - 1850
YEND_3000 = 4685 - 1850


def main():
    """Generate the experiment protocol figure."""

    # -------------------------------------------------------------------------
    # Load data
    # -------------------------------------------------------------------------

    ohc_pic = load_pic_ohc_1D()/ECHELLE_OHC
    
    amoc_pic = load_pic_amoc()

    ohc_hist_1000 = load_integrated_ohc("hist_tot", "0_btm")

    ohc_hist_3000 = load_integrated_ohc("hist_3000", "0_btm")

    starty_tot, starty_3000 = load_branching_years()

    # -------------------------------------------------------------------------
    # Figure setup
    # -------------------------------------------------------------------------
    cmap_1000 = get_cmap("Oranges")
    cmap_3000 = get_cmap("Blues")

    fig, ax1 = plt.subplots(figsize=(15, 10))

    # -------------------------------------------------------------------------
    # Reference OHC
    # -------------------------------------------------------------------------
    ax1.plot(
        ohc_pic.time,
        ohc_pic,
        color="black",
        linewidth=3,
        label="Reference simulation OHC",
    )

    # -------------------------------------------------------------------------
    # Historical members branched after 1000 years spin-up
    # -------------------------------------------------------------------------
    n_tot = len(starty_tot)

    for i, year in enumerate(starty_tot):
        color = cmap_1000(i / max(n_tot - 1, 1))

        ax1.plot(
            year - 1850,
            0,
            marker="o",
            color=color,
            markersize=10,
            transform=ax1.get_xaxis_transform(),
        )

        hist_member = ohc_hist_1000.isel(ensemble=i).assign_coords(time=TIME_HIST)

        ax1.plot(
            TIME_HIST + year - 1850,
            hist_member,
            "--",
            color=color,
            linewidth=2,
        )

    # -------------------------------------------------------------------------
    # Historical members branched after 3000 years spin-up
    # -------------------------------------------------------------------------
    n_3000 = len(starty_3000)

    for i, year in enumerate(starty_3000):
        color = cmap_3000(i / max(n_3000 - 1, 1))

        ax1.plot(
            year - 1850,
            0,
            marker="o",
            color=color,
            markersize=10,
            transform=ax1.get_xaxis_transform(),
        )

        hist_member = ohc_hist_3000.isel(ensemble=i).assign_coords(time=TIME_HIST)

        ax1.plot(
            TIME_HIST + year - 1850,
            hist_member,
            "--",
            color=color,
            linewidth=2,
        )

    # -------------------------------------------------------------------------
    # Shaded windows
    # -------------------------------------------------------------------------
    ax1.axvspan(
        YINI_1000,
        YEND_1000,
        color="orange",
        alpha=0.3,
        label="after 1000 years spin-up",
    )

    ax1.axvspan(
        YINI_3000,
        YEND_3000,
        color="blue",
        alpha=0.3,
        label="after 3000 years spin-up",
    )

    # -------------------------------------------------------------------------
    # Main axis formatting
    # -------------------------------------------------------------------------
    ax1.set_xlabel("Time (simulated years)", fontsize=30)
    ax1.set_ylabel("Global OHC (1e3 ZJ)", fontsize=30, color="black")
    ax1.tick_params(axis="x", labelsize=20, labelcolor="black")
    ax1.tick_params(axis="y", labelsize=20, labelcolor="black")
    ax1.grid(True)
    ax1.set_xlim(0, 3000)

    # -------------------------------------------------------------------------
    # Secondary axis: AMOC
    # -------------------------------------------------------------------------
    ax2 = ax1.twinx()
    ax2.plot(
        amoc_pic.year,
        amoc_pic,
        color="grey",
        linewidth=2,
        alpha=0.4,
        label="Reference simulation AMOC",
    )
    ax2.set_ylabel("AMOC (Sv)", fontsize=30, color="grey")
    ax2.tick_params(axis="y", labelsize=20, labelcolor="grey")

    ax2.set_zorder(0)
    ax1.set_zorder(1)
    ax1.patch.set_visible(False)

    # -------------------------------------------------------------------------
    # Custom legend
    # -------------------------------------------------------------------------
    legend_elements = [
        Patch(facecolor="orange", edgecolor="orange", label="hist+1000"),
        Patch(facecolor="blue", edgecolor="blue", label="hist+3000"),
        Line2D(
            [0], [0],
            marker="o",
            color="black",
            linestyle="None",
            markersize=10,
            label="Branching time",
        ),
        Line2D(
            [0], [0],
            color="black",
            lw=3,
            label="Reference simulation OHC",
        ),
        Line2D(
            [0], [0],
            color="grey",
            lw=2,
            alpha=0.4,
            label="Reference simulation AMOC",
        ),
    ]

    ax1.legend(
        handles=legend_elements,
        loc="center",
        bbox_to_anchor=(0.5, -0.2),
        ncol=2,
        fontsize=20,
    )

    ax1.set_title("Experiment protocol", fontsize=35)

    plt.tight_layout()

    # Save figure
    plt.savefig(FIG_DIR / "figure_1.pdf", bbox_inches="tight")


if __name__ == "__main__":
    main()
