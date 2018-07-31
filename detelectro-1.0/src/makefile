# detelectro -Determinant Quantum Monte Carlo for the Hubbard model
# Created by Francisco Brito July 25, 2018

#	The program is prepared to handle a variety of geometries (listed below).

# Input the number corresponding to the desired geometry
#		(1) 1d-pbc
#		(2) 1d-obc
#		(3) square-pbc
#		(4) square-obc
#		(5) rec-pbc
#		(6) rec-obc
#		(7) honeycomb-pbc
#		(8) hc-nanoribbon
#		(9)	hc-dot
# The NY parameter defines the width of the sample

# Set parameters of the simulation here.
CXX = g++ -DNSITES=4 -DDT_INV=8 -DBETA=6 -DHOPPING=1 -DGREEN_AFRESH_FREQ=4\
-DGEOM=1 -DNY=4
CXXFLAGS = -Wall -g -O2 -std=c++17 -I../includes

main: main.o
	$(CXX) $(CXXFLAGS) -o main main.o
	@echo "Arguments of main: U, mu, number of sweeps,\
	 number of warm-up sweeps, number of auto-correlation sweeps"

main.o: main.cpp ../includes/matrixgen.h ../includes/green.h ../includes/QDT.h\
	 ../includes/TDQ.h


clean:
	rm main main.o
