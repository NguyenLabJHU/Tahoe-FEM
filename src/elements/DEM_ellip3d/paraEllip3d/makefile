#########################################################
###        Part A: change with your environment       ###
#########################################################
###1. system platform
PLATFORM = supercomputer
#PLATFORM = workstation

###2. compilers package
#COMPILERS = gcc
COMPILERS = intel

###3. boost location
ifeq ($(PLATFORM), supercomputer)
 ifeq ($(COMPILERS), gcc)
  BOOST_ROOT = /work4/projects/openmpi/boost-1.55.0-openmpi-1.7.5-gcc-compilers-4.7.3
 endif
 ifeq ($(COMPILERS), intel)
  #BOOST_ROOT = /work4/projects/openmpi/boost-1.53.0-intelmpi-4.0.3.008
  BOOST_ROOT = /work4/projects/openmpi/boost-1.55.0-intelmpi-4.1.3.048-intel-compilers-13.0.1
 endif
endif

ifeq ($(PLATFORM), workstation)
 BOOST_ROOT = /usr/local/boost-1.53.0-openmpi-1.6.4-gcc-4.6.4
 #BOOST_ROOT = default
endif

#########################################################
###         Part B: usually do not need change        ###
#########################################################
###1. MPICXX, OPENMPI and OPTIMIZE
# C++ compiler
ifeq ($(COMPILERS), gcc)
 MPICXX = mpic++
endif

ifeq ($(COMPILERS), intel)
 MPICXX = mpiicpc
endif

# OpenMP
ifeq ($(COMPILERS), gcc)
 OPENMP = -fopenmp
endif

ifeq ($(COMPILERS), intel)
 OPENMP = -openmp -mt_mpi
endif

# optimize or debug
ifeq ($(COMPILERS), gcc)
 #full level of IEEE 754 compliance
 #OPTIMIZE = -O3 -DNDEBUG -frounding-math -fsignaling-nans
 OPTIMIZE = -O3 -DNDEBUG
endif
ifeq ($(COMPILERS), intel)
 #-fp-model source also implies keyword precise; and it is nearly as twice slow as default -fp-model fast=1
 #OPTIMIZE = -O3 -DNDEBUG -fp-model source
 OPTIMIZE = -O3 -DNDEBUG
endif
#OPTIMIZE = -Wall -g

# GCC quadmath, default off
#QUADMATH = -DQUADMATH -L/usr/local/gcc-4.6.2/lib64 -lquadmath

# CXXFLAGS
CXXFLAGS = $(OPTIMIZE) $(OPENMP) $(QUADMATH)

###2. BOOST details
ifeq ($(PLATFORM), supercomputer)
 BOOST_INCLUDE=-I$(BOOST_ROOT)/include
 BOOST_LIB_DIR=-L$(BOOST_ROOT)/lib
 BOOST_LIBS=-lboost_mpi -lboost_serialization -lboost_timer
 BOOST_RUN_LIB_DIR=-Wl,-rpath=$(BOOST_ROOT)/lib
endif

ifeq ($(PLATFORM), workstation)
  ifeq ($(BOOST_ROOT), default)
   BOOST_LIB_DIR=-L/usr/lib64/openmpi/lib:/usr/lib64/compat-openmpi/lib:/usr/lib64
   BOOST_LIBS=-lboost_mpi -lboost_serialization -lboost_timer
   BOOST_RUN_LIB_DIR=-Wl,-rpath=/usr/lib64/openmpi/lib:/usr/lib64/compat-openmpi/lib:/usr/lib64
  else  
   BOOST_INCLUDE=-I$(BOOST_ROOT)/include
   BOOST_LIB_DIR=-L$(BOOST_ROOT)/lib
   BOOST_LIBS=-lboost_mpi -lboost_serialization -lboost_timer
   BOOST_RUN_LIB_DIR=-Wl,-rpath=$(BOOST_ROOT)/lib
  endif
endif

INCLUDE=$(BOOST_INCLUDE)
LIB_DIR=$(BOOST_LIB_DIR)
LIBS=$(BOOST_LIBS)
RUN_LIB_DIR=$(BOOST_RUN_LIB_DIR)

###3. makefile
SOURCES = $(wildcard *.cpp)
OBJECTS = $(SOURCES:.cpp=.o)
#OBJECTS = $(patsubst %.cpp, %.o, $(SOURCES))
EXECUTABLE = paraEllip3d

.PHONY: all tar clean

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(MPICXX) -o $@ $(CXXFLAGS) $(OBJECTS) $(INCLUDE) $(LIB_DIR) $(LIBS) $(RUN_LIB_DIR)

%.o: %.cpp
	$(MPICXX) -c $< -o $@ $(CXXFLAGS) $(INCLUDE) $(LIB_DIR) $(LIBS)

sinclude $(SOURCES:.cpp=.d)

%.d: %.cpp
	$(MPICXX) -MM $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

tar:
	tar -cvf $(EXECUTABLE).tar *.h *.cpp makefile* readme

clean:
	-rm -f *.o *.d  *~ *.tar $(EXECUTABLE)

# sinclude is always resolved even if make tar/clean
