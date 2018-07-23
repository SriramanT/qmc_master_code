CXX = g++
CXXFLAGS = -Wall -g -std=c++11 -I.

hubbard-1d-finiteT: hubbard-1d-finiteT.o matrixgen.o prints.o green.o UDV.o VDU.o
	$(CXX) $(CXXFLAGS) -o hubbard-1d-finiteT hubbard-1d-finiteT.o matrixgen.o prints.o green.o UDV.o VDU.o

hubbard-1d-finiteT.o: hubbard-1d-finiteT.cpp matrixgen.cpp prints.cpp green.cpp matrixgen.h prints.h green.h

#hubbard-1d-finiteT-sign-problem: hubbard-1d-finiteT-sign-problem.o matrixgen.o prints.o green.o UDV.o VDU.o
#	$(CXX) $(CXXFLAGS) -o hubbard-1d-finiteT-sign-problem hubbard-1d-finiteT-sign-problem.o matrixgen.o prints.o green.o UDV.o VDU.o

#hubbard-1d-finiteT-sign-problem.o: hubbard-1d-finiteT-sign-problem.cpp matrixgen.cpp prints.cpp green.cpp matrixgen.h prints.h green.h
#	$(CXX) $(CXXFLAGS) -c hubbard-1d-finiteT-sign-problem.cpp

#tests: tests.o matrixgen.o prints.o green.o UDV.o VDU.o
#	$(CXX) $(CXXFLAGS) -o tests tests.o matrixgen.o prints.o green.o UDV.o VDU.o

#tests.o: tests.cpp matrixgen.cpp prints.cpp green.cpp VDU.cpp matrixgen.h prints.h green.h UDV.h VDU.h
#	$(CXX) $(CXXFLAGS) -c tests.cpp

matrixgen.o: matrixgen.h

prints.o: prints.h

UDV.o: UDV.h

VDU.o: VDU.h

green.o: green.cpp UDV.cpp VDU.cpp green.h UDV.h VDU.h
	$(CXX) $(CXXFLAGS) -c green.cpp
