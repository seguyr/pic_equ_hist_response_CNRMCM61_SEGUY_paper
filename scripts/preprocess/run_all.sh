#!/bin/bash


set -euo pipefail


SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"


echo "======================================"
echo "Running full preprocessing pipeline"
echo "======================================"


echo ""
echo "Step 1: OHC preprocessing"
echo "--------------------------------------"

bash ohc/run_all_members.sh 
bash ohc/04_build_hist_tot.sh 
bash ohc/05_build_hist_3000.sh 
bash ohc/06_compute_ohc_picontrol.sh 
bash ohc/07_build_pic_members.sh 
bash ohc/08_build_pic_tot.sh 
bash ohc/09_build_pic_3000.sh


echo ""
echo "Step 2: AMOC preprocessing"
echo "--------------------------------------"

bash amoc/11_build_pic_amoc_members.sh 
bash amoc/12_build_pic_amoc_tot.sh 
bash amoc/13_build_pic_amoc_3000.sh 
bash amoc/14_build_hist_amoc_members.sh 
bash amoc/15_build_hist_amoc_tot.sh 
bash amoc/16_build_hist_amoc_3000.sh 
bash amoc/17_build_picontrol_amoc_timeseries.sh 


echo ""
echo "======================================"
echo "Preprocessing completed successfully"
echo "======================================"
