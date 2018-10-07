#!/bin/bash

#PBS -l nodes=1:ppn=1
#PBS -l walltime=96:00:00

cd ~/qmc_master_code/DETELECTRO-1.0/fast-version/

make clean

make nsites=$1 beta=$2 source=$6

./simulation_$5_nsites$1_dt_inv16_beta$2_green_afresh_freq4_U_$3 1 $3 $4 15 $5 10448 256 2

cd results

echo "Saving simulation data"

python3 save-simulation-data.py
