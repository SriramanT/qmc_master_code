CXX = g++
CXXFLAGS = -Wall -g -std=c++11 -I.

main: main.o matrixgen.o prints.o green.o UDV.o VDU.o
	$(CXX) $(CXXFLAGS) -o main main.o matrixgen.o prints.o green.o UDV.o VDU.o

main.o: main.cpp matrixgen.cpp prints.cpp green.cpp matrixgen.h prints.h green.h

matrixgen.o: matrixgen.h

prints.o: prints.h

UDV.o: UDV.h

VDU.o: VDU.h

green.o: green.cpp UDV.cpp VDU.cpp green.h UDV.h VDU.h
	$(CXX) $(CXXFLAGS) -c green.cpp
