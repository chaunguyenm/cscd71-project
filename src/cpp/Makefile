CXX = g++
MPICXX = mpicxx
CXXFLAGS = -O2 -march=native -std=c++17 -g -Wall -Wunknown-pragmas
LDFLAGS = -g
OMPFLAGS = -fopenmp
GPROFFLAGS = -g -pg
INCLUDE = -I./include

all: main main-mpi

main: main.o fft.o
	$(CXX) $(INCLUDE) $(LDFLAGS) $(GPROFFLAGS) $(OMPFLAGS) -o $@ $^ $(LDLIBS)

main-mpi: main.o fft-mpi.o
	$(MPICXX) $(INCLUDE) $(LDFLAGS) $(GPROFFLAGS) $(OMPFLAGS) -o $@ $^ $(LDLIBS)

main.o: main.cc
	$(CXX) $(INCLUDE) $(GPROFFLAGS) $(OMPFLAGS) -c $(CXXFLAGS) -o $@ $<

fft.o: fft.cc fft.h
	$(CXX) $(INCLUDE) $(GPROFFLAGS) $(OMPFLAGS) -c $(CXXFLAGS) -o $@ $<

fft-mpi.o: fft.cc fft.h
	$(MPICXX) $(INCLUDE) $(GPROFFLAGS) $(OMPFLAGS) -c $(CXXFLAGS) -o $@ $<

clean:
	$(RM) main.o main-mpi.o fft.o fft-mpi.o

distclean: clean
	$(RM) main main-mpi

.PHONY: all clean distclean run 
