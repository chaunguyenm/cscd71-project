CXX = g++
MPICXX = mpicxx
CXXFLAGS = -O2 -march=native -std=c++17 -g -Wall -Wno-unknown-pragmas -Wno-sign-compare
LDFLAGS = -g
OMPFLAGS = -fopenmp
GPROFFLAGS = -g -pg
INCLUDE = -I./include

all: stft 

stft: stft.o fft.o
	$(MPICXX) $(INCLUDE) $(LDFLAGS) $(GPROFFLAGS) $(OMPFLAGS) -o $@ $^ $(LDLIBS)

stft.o: stft.cc
	$(MPICXX) $(INCLUDE) -c $(CXXFLAGS) $(GPROFLAGS) $(OMPFLAGS) -o $@ $<

fft.o: fft.cc fft.h
	$(MPICXX) $(INCLUDE) -c $(CXXFLAGS) $(GPROFFLAGS) $(OMPFLAGS) -o $@ $<

clean:
	$(RM) stft.o fft.o

distclean: clean
	$(RM) stft

.PHONY: all clean distclean 
