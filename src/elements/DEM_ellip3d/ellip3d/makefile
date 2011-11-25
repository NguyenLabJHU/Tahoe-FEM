# 1. You can toggle on/off OpenMP by uncomment/commenting the line
# 2. You can toggle on/off quadmath by uncomment/commenting the line

# C++ compiler
CXX = g++

# optimization
OPTIMIZE = -O3

# debugging options (macro DEBUG is defined in source code)
DEBUG = -Wall -g -DDEBUG

# OpenMP (macro OPENMP is defined in source code)
OPENMP = -DOPENMP -fopenmp

# quadmath
#QUADMATH = -DQUADMATH -L/usr/local/gcc-4.6.2/lib64 -lquadmath

# CXXFLAGS
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
