#!/bin/bash

#cd ~/qmc_master_code/DETELECTRO-1.0/fast-version/

make clean

make nsites=$1 beta=$2

for i in $(eval echo {0..$3})
  do
    ./simulation 1 $i 3.5 15 $4 10 5 1

    cd results

    echo "Saving simulation data"

    python save-simulation-data.py

    cd ../

  done
