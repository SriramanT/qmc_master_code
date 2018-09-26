#!/bin/bash

#cd ~/qmc_master_code/DETELECTRO-1.0/fast-version/

make clean

make nsites=$1 beta=$2

./simulation 1 $3 3.5 1 0 5096 512 8

cd results

echo "Saving simulation data"

python save-simulation-data.py
