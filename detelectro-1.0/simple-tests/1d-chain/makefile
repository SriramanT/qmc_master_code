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
# The NY parameter defines the width of the sample.
#	It is only meaningful for options 5, 6, 8, and 9.

nsites = 2
dt_inv = 64
beta = 2
hopping = 1
green_afresh_freq = 4
geom = 1
ny = 4

# Set parameters of the simulation here.
CXX = g++ -DNSITES=$(nsites) -DDT_INV=$(dt_inv) -DBETA=$(beta) -DHOPPING=$(hopping)\
 -DGREEN_AFRESH_FREQ=$(green_afresh_freq) -DGEOM=$(geom) -DNY=$(ny)

include_dir = ../../includes

CXXFLAGS = -Wall -g -O3 -std=c++11 -I$(include_dir)

main: main.o
	$(CXX) $(CXXFLAGS) -o main main.o
	@echo "Arguments of main: U, mu, number of sweeps"

main.o: main.cpp $(include_dir)/matrixgen.h $(include_dir)/green.h $(include_dir)/QDT.h\
	 $(include_dir)/TDQ.h

clean:
	rm -f main main.o
