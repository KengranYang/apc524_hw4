
objects =   heat_serial.o

CXXFLAGS = -g -Wall

all: heat_omp heat_serial

heat_omp: heat_omp.cc
	g++ -fopenmp -o $@ $^

heat_serial : $(objects)
	$(CXX) -o $@ $^



clean:
	$(RM) *.o
	$(RM) .depend

depend:
	$(CXX) -MM $(CXXFLAGS) *.cc > .depend

-include .depend
