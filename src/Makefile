# SYNOPSIS:
#
#   make [all]  	- makes everything.
#   make TARGET 	- makes the given target.
#   make clean  	- removes all files except executables.
#	make distclean 	- removes executables.

# Points to the root of Google Test, relative to where this file is.
GTEST_DIR=../googletest/googletest

# Where to find user code.
USER_DIR=.

# All Google Test headers.
GTEST_HEADERS = $(GTEST_DIR)/include/gtest/*.h \
                $(GTEST_DIR)/include/gtest/internal/*.h

# Compilers. Use MPICXX to compile MPI code and CXX otherwise.
CXX = g++
MPICXX = mpicxx

# Flags to compiler, defined as separate variables for ease of use.
CXXFLAGS += -O2 -march=native -std=c++17
CXXFLAGS += -g -Wall -Wextra -Wno-unknown-pragmas -Wno-sign-compare
OMPFLAGS = -fopenmp
GPROFFLAGS = -g -pg
FFTWFLAGS = -I/usr/local/include -L/usr/local/lib -lfftw3 -lm

# Flags passed to the preprocessor.
# Set Google Test's header directory as a system directory, such that
# the compiler doesn't generate warnings in Google Test headers.
CPPFLAGS += -isystem $(GTEST_DIR)/include

# Path for linking.
INCLUDE = -I./include

# House-keeping build targets.
all: stft 

clean:
	$(RM) *.o gtest.a gtest_main.a

distclean:
	$(RM) fft_unittest stft

# Builds gtest.a and gtest_main.a.
# Usually you shouldn't tweak such internal variables, indicated by a
# trailing _.
GTEST_SRCS_ = $(GTEST_DIR)/src/*.cc $(GTEST_DIR)/src/*.h $(GTEST_HEADERS)

# For simplicity and to avoid depending on Google Test's
# implementation details, the dependencies specified below are
# conservative and not optimized.  This is fine as Google Test
# compiles fast and for ordinary users its source rarely changes.
gtest-all.o : $(GTEST_SRCS_)
	$(CXX) $(CPPFLAGS) -I$(GTEST_DIR) $(CXXFLAGS) -c \
            $(GTEST_DIR)/src/gtest-all.cc

gtest_main.o : $(GTEST_SRCS_)
	$(CXX) $(CPPFLAGS) -I$(GTEST_DIR) $(CXXFLAGS) -c \
            $(GTEST_DIR)/src/gtest_main.cc

gtest.a : gtest-all.o
	$(AR) $(ARFLAGS) $@ $^

gtest_main.a : gtest-all.o gtest_main.o
	$(AR) $(ARFLAGS) $@ $^	

# Builds the project.
stft: stft.o fft.o
	$(MPICXX) $(INCLUDE) $(GPROFFLAGS) $(OMPFLAGS) -o $@ $^ $(LDLIBS)

stft.o: stft.cc
	$(MPICXX) $(INCLUDE) -c $(CXXFLAGS) $(GPROFLAGS) $(OMPFLAGS) -o $@ $<

fft.o: fft.cc fft.h
	$(MPICXX) $(INCLUDE) -c $(CXXFLAGS) $(GPROFFLAGS) $(OMPFLAGS) -o $@ $<		

# Builds unit test.  A test should link with either gtest.a or
# gtest_main.a, depending on whether it defines its own main()
# function.
fft_unittest.o : $(USER_DIR)/fft_unittest.cc $(GTEST_HEADERS)
	$(CXX) $(INCLUDE) $(CPPFLAGS) $(CXXFLAGS) $(FFTWFLAGS) -c $(USER_DIR)/fft_unittest.cc

fft_unittest : fft_unittest.o fft.o gtest_main.a
	$(MPICXX) $^ $(FFTWFLAGS) $(INCLUDE) $(CPPFLAGS) $(CXXFLAGS) $(OMPFLAGS) -o $@	

test: fft_unittest
	./fft_unittest

.PHONY: all clean distclean test
