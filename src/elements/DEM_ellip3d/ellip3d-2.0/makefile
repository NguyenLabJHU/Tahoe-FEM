# options available in this makefile:
# 1. option to toggle between gcc and intel
# 2. option to toggle on/off OpenMP (default on)
# 3. option to switch between optimization/debug mode (default optimization)
# 4. option to toggle on/off GCC quadmath (default off)

#1. gcc or intel
#COMPILERS = gcc
COMPILERS = intel

# C++ compiler
ifeq ($(COMPILERS), gcc)
 CXX = g++
endif

ifeq ($(COMPILERS), intel)
 CXX = icpc
endif

#2. OpenMP or not (macro OPENMP is defined in source code)
OMP=-DOPENMP
#OMP=

ifdef OMP
 ifeq ($(COMPILERS), gcc)
  OPENMP = $(OMP) -fopenmp
 endif
 ifeq ($(COMPILERS), intel)
  OPENMP = $(OMP) -openmp
 endif
endif

#3. optimize or debug
ifeq ($(COMPILERS), gcc)
 OPTIMIZE = -O3 -DNDEBUG
endif
ifeq ($(COMPILERS), intel)
 #-fp-model source also implies keyword precise; and it is much slower than default -fp-model fast=1
 #OPTIMIZE = -O3 -DNDEBUG -fp-model source
 OPTIMIZE = -O3 -DNDEBUG
endif
#OPTIMIZE = -Wall -g

#4. GCC quadmath
#QUADMATH = -DQUADMATH -L/usr/local/gcc-4.6.2/lib64 -lquadmath

CXXFLAGS = $(OPTIMIZE) $(OPENMP) $(QUADMATH)

SOURCES = $(wildcard *.cpp)
OBJECTS = $(SOURCES:.cpp=.o)
#OBJECTS = $(patsubst %.cpp, %.o, $(SOURCES))
EXECUTABLE = ellip3d

.PHONY: all tar clean

# template classes contain implementation and need to recompile upon change
all: $(EXECUTABLE) contact.h boundary.h

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $@

%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) $< -o $@

sinclude $(SOURCES:.cpp=.d)

%.d: %.cpp
	$(CXX) -MM $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

tar:
	tar -cvf $(EXECUTABLE).tar *.h *.cpp makefile* readme
clean:
	-rm -f *.o *.d  *~ *.tar $(EXECUTABLE)

# sinclude is always resolved even if make tar/clean
