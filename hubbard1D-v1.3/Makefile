CXX = g++
CXXFLAGS = -Wall -g -std=c++11 -I.

hubbard-1d-finiteT: hubbard-1d-finiteT.o matrixgen.o prints.o green.o UDV.o
	$(CXX) $(CXXFLAGS) -o hubbard-1d-finiteT hubbard-1d-finiteT.o matrixgen.o prints.o green.o UDV.o

hubbard-1d-finiteT.o: hubbard-1d-finiteT.cpp matrixgen.cpp prints.cpp green.cpp matrixgen.h prints.h green.h
	$(CXX) $(CXXFLAGS) -c hubbard-1d-finiteT.cpp

matrixgen.o: matrixgen.h

prints.o: prints.h

UDV.o: UDV.h

green.o: UDV.cpp green.h UDV.h
	$(CXX) $(CXXFLAGS) -c green.cpp
