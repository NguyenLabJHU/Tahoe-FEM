# options available in this makefile:
# 1. option to toggle on/off OpenMP (default on)
# 2. option to toggle on/off GCC quadmath (default off)
# 3. option to switch between optimization/debug mode (default optimization)

# C++ compiler
CXX = g++

# optimization
OPTIMIZE = -O3

# debugging options (macro DEBUG is defined in source code)
DEBUG = -Wall -g -DDEBUG

# 1. OpenMP (macro OPENMP is defined in source code)
OPENMP = -DOPENMP -fopenmp

# 2. GCC quadmath
#QUADMATH = -DQUADMATH -L/usr/local/gcc-4.6.2/lib64 -lquadmath

# 3. Optimize or debug
CXXFLAGS = $(OPTIMIZE) $(OPENMP) $(QUADMATH)
#CXXFLAGS = $(DEBUG) $(OPENMP) $(QUADMATH)

SOURCES = $(wildcard *.cpp)
OBJECTS = $(SOURCES:.cpp=.o)
#OBJECTS = $(patsubst %.cpp, %.o, $(SOURCES))
EXECUTABLE = ellip3d

.PHONY: all clean

all: $(SOURCES) $(EXECUTABLE)

sinclude $(SOURCES:.cpp=.d)

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $@

%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) $< -o $@

%.d: %.cpp
	$(CXX) -MM $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

clean:
	-rm -f *.o *.d *~ ellip3d
