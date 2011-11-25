# 1. You can turn on/off openmp

# C++ compiler
CXX = g++

# optimization
OPTIMIZE = -O3

# debugging options (macro DEBUG is defined in source code)
DEBUG = -Wall -g -DDEBUG

# OPENMP (macro OPENMP is defined in source code)
OPENMP = -DOPENMP -fopenmp
#OPENMP = 

# CXXFLAGS
 CXXFLAGS = $(OPTIMIZE) $(OPENMP)
#CXXFLAGS = $(DEBUG) $(OPENMP)

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
