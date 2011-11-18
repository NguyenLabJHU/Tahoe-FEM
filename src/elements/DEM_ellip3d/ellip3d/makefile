# This is the makefile for ellip3d. There are 4 modes to select from:
# 1. optimization                                                    
# 2. optimization + openmp                                           
# 3. debug                                                           
# 4. debug + openmp                                                  
# What you need to do is uncommenting one of the 4 lines enclosed 
# below, and you are only recommended to select mode 1 or 2. 
# Mode 3 and 4 write massive debug information on disk so they are 
# only applicable to cases of a few particles or a few time steps. 

# C++ compiler
CXX = g++

# optimization
OPTIMIZE = -O3

# debugging options (macro DEBUG is defined in source code)
DEBUG = -Wall -g -DDEBUG

# OPENMP (macro OPENMP is defined in source code)
OPENMP = -DOPENMP -fopenmp

######################################################################
#CXXFLAGS = $(OPTIMIZE)
 CXXFLAGS = $(OPTIMIZE) $(OPENMP)
#CXXFLAGS = $(DEBUG)
#CXXFLAGS = $(DEBUG) $(OPENMP)
######################################################################

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
