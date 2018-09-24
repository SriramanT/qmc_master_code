#!/bin/bash

#cd ~/qmc_master_code/DETELECTRO-1.0/fast-version/

make clean

make nsites=$1 beta=$2

for i in $(eval echo {1..$4})
  do
    ./simulation 1 $3 3.5 15 $((2**$i)) 10 5 1

    cd results

    echo "Saving simulation data"

    python save-simulation-data.py

    cd ../

  done
