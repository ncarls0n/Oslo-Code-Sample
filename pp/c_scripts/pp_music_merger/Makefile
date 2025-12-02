FC = gfortran
CXX = g++
FFLAGS = -fPIC
CXXFLAGS = -fPIC

all: main

array_module.o: array_module.f90
	$(FC) $(FFLAGS) -c array_module.f90

main.o: main.cpp
	$(CXX) $(CXXFLAGS) -c main.cpp

main: array_module.o main.o
	$(CXX) -o main array_module.o main.o

clean:
	rm -f *.mod *.o main
