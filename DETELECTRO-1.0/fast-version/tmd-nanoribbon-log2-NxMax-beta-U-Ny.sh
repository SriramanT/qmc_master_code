#!/bin/bash

#cd ~/qmc_master_code/DETELECTRO-1.0/fast-version/

for i in $(eval echo {1..$1})
  do
    make clean

    make nsites=$((3*$4*((2**$i)))) beta=$2

    ./simulation 1 $3 3.5 15 $4 10 5 1

    cd results

    echo "Saving simulation data"

    python save-simulation-data.py

    cd ../

  done
