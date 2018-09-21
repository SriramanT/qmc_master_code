#!/bin/bash

cd ~/qmc_master_code/DETELECTRO-1.0/fast-version/

make clean

make nsites=48 beta=16

./simulation 1 1 12.6 15 4 100 20 2

cd results

echo "Saving simulation data"

python save-simulation-data.py
