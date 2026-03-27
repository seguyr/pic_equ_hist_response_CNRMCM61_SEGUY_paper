"""
Figure 5: Historical time series response comaprison
=============================

This script plots:
- pic_tot time series in anomalies to 1850 - 1900, 
- pic_3000 time series in anomalies to 1850 - 1900, 
- hist_tot time series in anomalies to 1850 - 1900, 
- hist_3000 time series in anomalies to 1850 - 1900,
- the difference between historical ohc response dedrift by time matching

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
    time_matching,
    anomalies,
    gain,
    bootstrap,
    bootstrap_2,
    boot_1D, 
    boot_2D,
    get_stats,
    load_pic_ohc,
    load_pic_amoc,
    load_integrated_ohc,
    load_branching_years,
)

# -----------------------------------------------------------------------------
# Physical constants
# -----------------------------------------------------------------------------
ECHELLE_OHC = 1e21  # ZJ
