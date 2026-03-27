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

# -----------------------------------------------------------------------------
# Calcul
# -----------------------------------------------------------------------------

# Load 2D OHC ensembles
OHC_2D_hist_tot, OHC_2D_hist_3000, OHC_2D_pic_tot, OHC_2D_pic_3000 = load_ohc_2d_ensembles()

# dedrift time matching
ohc_dd_1000 = time_matching(OHC_2D_hist_tot, OHC_2D_pic_tot)
ohc_dd_3000 = time_matching(OHC_2D_hist_3000, OHC_2D_pic_3000)

# Anomalies relative to first 50 years
OHC_dd_1000 = anomalies(ohc_dd_1000)
OHC_dd_3000 = anomalies(ohc_dd_3000)

# Mean gain over final 20 years
m_OHC_dd_1000 = gain(OHC_dd_1000)
m_OHC_dd_3000 = gain(OHC_dd_3000)

# Bootstrap 
hist_cor_anom_1000 = boot(OHC_dd_1000, n_boot=n_boot)
hist_cor_anom_3000 = boot(OHC_dd_3000, n_boot=n_boot)
diff_hist_cor = boot_diff(OHC_dd_1000, OHC_dd_3000, n_boot=n_boot)
m_hist_cor_1000 = boot(m_OHC_dd_1000, n_boot=n_boot)
m_hist_cor_3000 = boot(m_OHC_dd_3000, n_boot=n_boot)
m_diff_hist_cor = boot_diff(m_OHC_dd_1000, m_OHC_dd_3000, n_boot=n_boot)






