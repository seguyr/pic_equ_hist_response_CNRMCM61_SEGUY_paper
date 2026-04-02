# pic_equ_hist_response_CNRMCM61_SEGUY_paper
# Reproducing the Figures from the Article

This repository contains the scripts and processed datasets required to reproduce the figures presented in the associated article.

### Repository Structure Overview

The workflow follows three main steps:

1. Download raw datasets
2. Run preprocessing scripts
3. Generate figures

---

## Generating the Figures

To reproduce the figures from the paper:

Go to the directory: /scripts

Run the Python script corresponding to the figure you want to generate.

Each script automatically creates the associated figures in: figures/.pdf  and figures/.png

---

## Data Used for Plotting

The final datasets used directly by the plotting scripts are located in: data/data_plot/

These correspond to the last processing step before figure generation.

---

## Metadata

Metadata required for plotting (masks, grids, branching years, etc.) are stored in: data/metadata/

---

## Raw Data Availability

Raw datasets are too large to be stored in this repository.

They are available here: https://esgf.github.io/index.html

This archive includes:

- CMIP6 ensemble members
- the 2000-year piControl simulations

Additional datasets used in the study:

- extended piControl simulations
- reconstructed ensembles

are available upon request from the authors.

---

## Preprocessing Workflow

After downloading the raw datasets, preprocessing must be performed before figure generation.
All preprocessing scripts are located in: scripts/preprocess/

To reproduce the processed datasets automatically, run: run_all.sh
This script executes all preprocessing steps in the correct order and recreates the intermediate datasets in the specified output directory.

It reproduces datasets corresponding to those partially available in:
data/
├── pic_tot
├── hist_tot
├── pic_3000
├── hist_3000
└── pic_global

---

## Intermediate Processing Functions

All functions used to transform intermediate datasets into the final plotting datasets are located in: functions/utils.py

These functions are used to generate: data/data_plot/

---

## Python Environment

The Python environment used for this project is provided here: pip install -r environment/requirements.txt

---






