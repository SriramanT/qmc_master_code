#!/bin/bash

#cd ~/qmc_master_code/ANALYTICAL/src

for i in {5..10}
  do

    python3 tmdnr-mf-Nx-Ny-beta-U-init.py $((2**$i)) 32 20 8 1

  done
