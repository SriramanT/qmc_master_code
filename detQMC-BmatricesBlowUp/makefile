CXX = g++
CXXFLAGS = -Wall -g -std=c++11 -I.

main: main.o matrixgen.o green.o UDV.o VDU.o
	$(CXX) $(CXXFLAGS) -o main main.o matrixgen.o green.o UDV.o VDU.o

main.o: main.cpp matrixgen.h green.h

matrixgen.o: matrixgen.h

UDV.o: UDV.h

VDU.o: VDU.h

green.o: green.h
