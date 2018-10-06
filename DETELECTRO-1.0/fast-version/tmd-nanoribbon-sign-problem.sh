#!/bin/bash

#cd ~/qmc_master_code/DETELECTRO-1.0/fast-version/

make clean

make nsites=$1 beta=$2 eq_or_uneq=src/mainSignProblem.cpp object=src/mainSignProblem.o

for i in $(eval echo {1..$4})
  do

    ./simulation 1 $3 $(echo 0.1*$i - 0.1 | bc) 15 $5 $6 $7 $8

    cd results

    echo "Saving simulation data"

    python3 save-sign.py

    cd ../

  done
