
objects =   heat_serial.o

CXXFLAGS = -g -Wall -std=c++0x
OMPFLAGS = -std=c++0x -fopenmp

all: heat_omp heat_serial heat_mpi

heat_omp: heat_omp.cc
	g++ $(OMPFLAGS) -o $@ $^

heat_serial : $(objects)
	$(CXX) -o $@ $^

heat_mpi:heat_mpi.cc
	mpicc -o $@ $^

clean:
	$(RM) *.o
	$(RM) .depend

depend:
	$(CXX) -MM $(CXXFLAGS) *.cc > .depend

-include .depend
