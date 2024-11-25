CXX = g++
CXXFLAGS = -O2 -march=native -std=c++17 -g -Wall
LDFLAGS = -g
OMPFLAGS = -fopenmp
GPROFFLAGS = -g -pg
OBJS = main.o fft.o
INCLUDE = -I./include

all: main 

main.o: main.cc
	$(CXX) $(INCLUDE) -c $(CXXFLAGS) -o $@ $<

fft.o: fft.cc fft.h
	$(CXX) $(INCLUDE) -c $(CXXFLAGS) -o $@ $<

main: ${OBJS}
	$(CXX) $(INCLUDE) $(LDFLAGS) -o $@ $^ $(LDLIBS)

clean:
	$(RM) $(OBJS)

distclean: clean
	$(RM) traffic traffic-????.npy traffic.png test-????.npy largetest-????.npy test.png largetest.png

.PHONY: all clean distclean run 
