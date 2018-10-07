# DETELECTRO - Determinant Quantum Monte Carlo for Hubbard-like models
# Created by Francisco Brito August 14, 2018

#	The program is prepared to handle a variety of geometries (listed below).

#		Input the number corresponding to the desired geometry
#		(1)		1D Periodic Chain
#		(2) 	1D Open Chain
#		(3) 	2D Periodic Square Lattice
#		(4) 	2D Open Square Lattice
#		(5) 	2D Periodic Rectangular Lattice
#		(6) 	2D Open Rectangular Lattice
#		(7) 	2D Periodic Triangular Lattice
#		(8) 	2D Nanoribbon Triangular Lattice
#		(9) 	2D Periodic Honeycomb Lattice
#		(10)	2D Honeycomb Nanoribbon
#		(11)	2D Honeycomb Strained Nanoribbon
#		(12)	2D Honeycomb Triangular Dot
#		(13)	2D Honeycomb Rectangular Dot
#		(14)	2D Minimal model of a periodic TMD sample (Liu2013)
#		(15)	2D Minimal model of a TMD nanoribbon (Liu2013)

#	DEFAULT PARAMETERS

nsites=2
# orbital + real space
dt_inv=16
beta=2
green_afresh_freq=4
U=20

source=mainEqTime

# eq_or_uneq=src/mainUneqTime.cpp
# object=src/mainUneqTime.o

# Set parameters of the simulation here.
CXX = g++ -DNSITES=$(nsites) -DDT_INV=$(dt_inv) -DBETA=$(beta)\
 -DGREEN_AFRESH_FREQ=$(green_afresh_freq)

include_dir = ../includes

CXXFLAGS = -Wall -g -O3 -std=c++11 -I$(include_dir)

simulation_$(source)_nsites$(nsites)_dt_inv$(dt_inv)_beta$(beta)_green_afresh_freq$(green_afresh_freq)_U_$(U): src/$(source).o
	$(CXX) $(CXXFLAGS) -o simulation_$(source)_nsites$(nsites)_dt_inv$(dt_inv)_beta$(beta)_green_afresh_freq$(green_afresh_freq)_U_$(U) src/$(source).o
	@echo ""
	@echo "		DETELECTRO - Determinant Quantum Monte Carlo for Hubbard-like models"
	@echo ""
	@echo "			Created by Francisco Brito (2018)"
	@echo ""
	@echo "The code has compiled successfully. To change the number of sites, \
	inverse Trotter error, inverse temperature, or the frequency of recomputing \
	the Green's functions, type:"
	@echo ""
	@echo "make clean"
	@echo ""
	@echo "make nsites=<Number of sites> dt_inv=<Inverse Trotter Error> \
	beta=<Inverse Temperature> green_afresh_freq=<Frequency of Recomputing G>"
	@echo ""
	@echo "To run a simulation, simply type ./simulation followed by its arguments:"
	@echo ""
	@echo "./simulation <t> <U> <mu> <geom> <Ny>\
	 <Total Number of Sweeps (Space-Time)>\
	 <Number of Warm-up Sweeps (Space-Time)> \
	 <Number of Auto-correlation Sweeps (Space-Time)>"
	@echo ""
	@echo "The program is prepared to handle a variety of geometries (listed below).\
	 Input the number corresponding to the desired geometry:"
	@echo ""
	@echo "(1)		1D Periodic Chain"
	@echo ""
	@echo "(2) 		1D Open Chain"
	@echo ""
	@echo "(3) 		2D Periodic Square Lattice"
	@echo ""
	@echo "(4) 		2D Open Square Lattice"
	@echo ""
	@echo "(5) 		2D Periodic Rectangular Lattice"
	@echo ""
	@echo "(6) 		2D Open Rectangular Lattice"
	@echo ""
	@echo "(7) 		2D Periodic Triangular Lattice"
	@echo ""
	@echo "(8) 		2D Nanoribbon Triangular Lattice"
	@echo ""
	@echo "(9) 		2D Periodic Honeycomb Lattice"
	@echo ""
	@echo "(10)		2D Honeycomb Nanoribbon"
	@echo ""
	@echo "(11)		2D Honeycomb Hexagonal Dot"
	@echo ""
	@echo "(12)		2D Honeycomb Triangular Dot"
	@echo ""
	@echo "(13)		2D Honeycomb Rectangular Dot"
	@echo ""
	@echo "(14)		2D Minimal model of a periodic TMD sample (Liu et al., Phys Rev B 88, 085433, 2013 )\
	 - nsites includes orbital space, i.e. nsites=n_orbitals * n_spatial_sites."
	@echo ""
	@echo "(15)		2D Minimal model of a TMD nanoribbon (Liu et al., Phys Rev B 88, 085433, 2013 )\
	 - nsites includes orbital space, i.e. nsites=n_orbitals * n_spatial_sites."
	@echo ""
	@echo "The Ny parameter is only meaningful for geometry options 5, 6, 8, 10, 13, and 15"
	@echo ""

src/$(source).o: src/$(source).cpp $(include_dir)/matrixgen.h $(include_dir)/green.h $(include_dir)/QDT.h\
	 $(include_dir)/TDQ.h

clean:
	rm -f simulation_$(source)_nsites$(nsites)_dt_inv$(dt_inv)_beta$(beta)_green_afresh_freq$(green_afresh_freq)_$(U) src/$(source).o
	rm -f ./temp-data/*
