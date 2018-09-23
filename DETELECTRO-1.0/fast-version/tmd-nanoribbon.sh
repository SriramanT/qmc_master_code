#!/bin/bash

#cd ~/qmc_master_code/DETELECTRO-1.0/fast-version/

make clean

make nsites=48 beta=16

./simulation 1 4 3.5 15 4 2048 512 2

cd results

echo "Saving simulation data"

python save-simulation-data.py
