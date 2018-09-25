#!/bin/bash

#cd ~/qmc_master_code/ANALYTICAL/src

for i in {3..7}
  do

    python3 tmdnr-mf-Nx-Ny-beta-U-init.py 512 $((2**$i)) 20 8 1

  done
