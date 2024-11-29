CXX = g++
MPICXX = mpicxx
CXXFLAGS = -O2 -march=native -std=c++17 -g -Wall -Wno-unknown-pragmas
LDFLAGS = -g
OMPFLAGS = -fopenmp
GPROFFLAGS = -g -pg
INCLUDE = -I./include

all: stft stft-omp

stft: stft.o fft.o
	$(MPICXX) $(INCLUDE) $(LDFLAGS) -o $@ $^ $(LDLIBS)

stft-omp: stft.o fft-omp.o
	$(MPICXX) $(INCLUDE) $(LDFLAGS) $(OMPFLAGS) -o $@ $^ $(LDLIBS)

stft.o: stft.cc
	$(MPICXX) $(INCLUDE) -c $(CXXFLAGS) -o $@ $<

fft.o: fft.cc fft.h
	$(MPICXX) $(INCLUDE) -c $(CXXFLAGS) -o $@ $<

fft-omp.o: fft.cc fft.h
	$(MPICXX) $(INCLUDE) $(OMPFLAGS) -c $(CXXFLAGS) -o $@ $<

clean:
	$(RM) stft.o fft.o fft-omp.o

distclean: clean
	$(RM) stft stft-omp

.PHONY: all clean distclean 
