#!/bin/bash

#PBS -l nodes=1:ppn=1
#PBS -l walltime=96:00:00

cd ~/qmc_master_code/DETELECTRO-1.0/fast-version/

make clean

make nsites=$1 beta=$2

./simulation 1 $3 3.5 15 $4 2048 512 4

cd results

echo "Saving simulation data"

python save-simulation-data.py
