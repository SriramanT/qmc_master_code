#!/bin/bash

cd ~/qmc_master_code/DETELECTRO-1.0/fast-version

make clean

make nsites=96 beta=8

./simulation 1 8 12.55 15 4 5012 512 4

cd results

echo "Saving simulation data"

python save-simulation-data.py
