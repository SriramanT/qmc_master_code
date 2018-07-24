CXX = g++
CXXFLAGS = -Wall -g -O2 -std=c++17 -I../includes

main: main.o
	$(CXX) $(CXXFLAGS) -o main main.o

main.o: main.cpp ../includes/matrixgen.h ../includes/green.h ../includes/QDT.h ../includes/TDQ.h
