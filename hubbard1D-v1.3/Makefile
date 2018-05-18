CXX = g++
CXXFLAGS = -Wall -g -std=c++11 -I.

hubbard-1d-finiteT: hubbard-1d-finiteT.o matrixgen.o prints.o green.o UDV.o
	$(CXX) $(CXXFLAGS) -o hubbard-1d-finiteT hubbard-1d-finiteT.o matrixgen.o prints.o green.o UDV.o

hubbard-1d-finiteT.o: hubbard-1d-finiteT.cpp matrixgen.cpp prints.cpp green.cpp matrixgen.h prints.h green.h
	$(CXX) $(CXXFLAGS) -c hubbard-1d-finiteT.cpp

#hubbard-1d-finiteT-sign-problem: hubbard-1d-finiteT-sign-problem.o matrixgen.o prints.o green.o UDV.o

#hubbard-1d-finiteT-sign-problem.o: hubbard-1d-finiteT-sign-problem.cpp matrixgen.cpp prints.cpp green.cpp matrixgen.h prints.h green.h
	$(CXX) $(CXXFLAGS) -c hubbard-1d-finiteT-sign-problem.cpp

#stableGreen: stableGreen.o matrixgen.o prints.o green.o UDV.o
#	$(CXX) $(CXXFLAGS) -o stableGreen stableGreen.o matrixgen.o prints.o green.o UDV.o

#stableGreen.o: stableGreen.cpp matrixgen.cpp prints.cpp green.cpp matrixgen.h prints.h green.h
#	$(CXX) $(CXXFLAGS) -c stableGreen.cpp

matrixgen.o: matrixgen.h

prints.o: prints.h

UDV.o: UDV.h

green.o: green.cpp UDV.cpp green.h UDV.h
#	$(CXX) $(CXXFLAGS) -c green.cpp
