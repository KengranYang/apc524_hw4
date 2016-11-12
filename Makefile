
objects =   hw4_serial.o

CXXFLAGS = -g -Wall

all: hw4_serial

hw4_serial : $(objects)
	$(CXX) -o $@ $^

clean:
	$(RM) *.o
	$(RM) .depend

depend:
	$(CXX) -MM $(CXXFLAGS) *.cc > .depend

-include .depend
