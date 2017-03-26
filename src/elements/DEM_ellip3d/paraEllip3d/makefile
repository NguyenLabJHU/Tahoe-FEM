#########################################################
###        Part A: change with your environment       ###
#########################################################
###1. system platform
PLATFORM = soilblast
#PLATFORM = spirit
#PLATFORM = garnet
#PLATFORM = janus
#PLATFORM = lightning
#PLATFORM = excalibur
#PLATFORM = thunder
#PLATFORM = topaz

###2. compiler combinations (c++ and mpi)
COMPILERS = gcc-openmpi
#COMPILERS = intel-intel
#COMPILERS = intel-sgi
#COMPILERS = intel-cray
#COMPILERS = gcc-sgi
#COMPILERS = gcc-cray

###3. compile papi for performance testing?
PAPI = no
#PAPI = yes

###4. boost location
ifeq ($(PLATFORM), soilblast)
 BOOST_ROOT = /usr/local/boost-1.53.0-openmpi-1.6.4-gcc-4.6.4
 #BOOST_ROOT = default
 ifeq ($(PAPI), yes)
  PAPI_ROOT = /usr/local/papi-5.4.1-openmpi-1.6.4-gcc-4.6.4
 endif
endif

ifeq ($(PLATFORM), spirit)
 ifeq ($(COMPILERS), gcc-openmpi)
  BOOST_ROOT = /work4/projects/openmpi/boost-1.55.0-openmpi-1.7.5-gcc-compilers-4.7.3
 endif
 ifeq ($(COMPILERS), intel-intel)
  BOOST_ROOT = /work4/projects/openmpi/boost-1.55.0-intelmpi-4.1.3.048-intel-compilers-13.0.1
 endif
 ifeq ($(COMPILERS), intel-sgi)
  BOOST_ROOT = /home/yanb/local/boost-1.57.0_mpt-2.14_intel-compilers-16.0.3
  ifeq ($(PAPI), yes)
   PAPI_ROOT = /home/yanb/local/papi-5.4.1_mpt-2.14_intel-compilers-16.0.3
  endif
 endif
endif

ifeq ($(PLATFORM), garnet)
 ifeq ($(COMPILERS), intel-cray)
  BOOST_ROOT = /u/yanb/local/boost-1.55.0_cray-mpich6.0.0_PrgEnv-intel4.2.34
 endif
endif

ifeq ($(PLATFORM), lightning)
 ifeq ($(COMPILERS), intel-cray)
  BOOST_ROOT = /home/yanb/local/boost-1.55.0_cray-mpich7.1.3_PrgEnv-intel5.2.40
 endif
endif

ifeq ($(PLATFORM), excalibur)
 ifeq ($(COMPILERS), intel-cray)
  BOOST_ROOT = /usr/people/yanb/local/boost-1.57.0_cray-mpich7.1.0_PrgEnv-intel5.2.40
  ifeq ($(PAPI), yes)
   PAPI_ROOT = /usr/people/yanb/local/papi-5.4.1_cray-mpich7.1.0_PrgEnv-intel5.2.40
  endif
 endif
 ifeq ($(COMPILERS), gcc-cray)
  BOOST_ROOT = /usr/people/yanb/local/boost-1.57.0_cray-mpich7.1.0_PrgEnv-gnu5.2.40
 endif
 ifeq ($(COMPILERS), gcc-openmpi)
  BOOST_ROOT = /usr/people/yanb/local/boost-1.63.0_openmpi-1.8.4_PrgEnv-gnu5.2.40
  ifeq ($(PAPI), yes)
   PAPI_ROOT = /usr/people/yanb/local/papi-5.5.1_openmpi-1.8.4_PrgEnv-gnu5.2.40 
  endif
 endif
endif

ifeq ($(PLATFORM), thunder)
 ifeq ($(COMPILERS), intel-sgi)
  BOOST_ROOT = /home/yanb/local/boost-1.57.0_mpt-2.12_intel-compilers-15.3.187
  ifeq ($(PAPI), yes)
   PAPI_ROOT = /home/yanb/local/papi-5.4.1_mpt-2.12_intel-compilers-15.3.187
  endif
 endif
 ifeq ($(COMPILERS), gcc-sgi)
  BOOST_ROOT = /home/yanb/local/boost-1.57.0_mpt-2.12_gcc-compilers-4.4.6
 endif
endif

ifeq ($(PLATFORM), topaz)
 ifeq ($(COMPILERS), intel-sgi)
  BOOST_ROOT = /p/home/yanb/local/boost-1.57.0_mpt-2.14_intel-16.0.0
  ifeq ($(PAPI), yes)
   PAPI_ROOT = /p/home/yanb/local/papi-5.4.1_mpt-2.14_intel-16.0.0
  endif
 endif
endif

ifeq ($(PLATFORM), janus)
 ifeq ($(COMPILERS), gcc-openmpi)
  BOOST_ROOT = /home/yanb/local/boost-1.55.0_openmpi-1.8.2_gcc-4.9.1
 endif
 ifeq ($(COMPILERS), intel-intel)
  BOOST_ROOT = /home/yanb/local/boost-1.55.0_impi-13.0.0_intel-13.0.0
 endif
endif

#########################################################
###         Part B: usually do not change             ###
#########################################################
###1. MPI wrapper, OpenMP and optimize/debug
# MPI wrapper
ifeq ($(COMPILERS), gcc-openmpi)
 MPICXX = mpic++
endif

ifeq ($(COMPILERS), gcc-sgi)
 MPICXX = mpicxx
endif

ifeq ($(COMPILERS), gcc-cray)
 MPICXX = CC
endif

ifeq ($(COMPILERS), intel-intel)
 MPICXX = mpiicpc
endif

ifeq ($(COMPILERS), intel-sgi)
 MPICXX = mpicxx
endif

ifeq ($(COMPILERS), intel-cray)
 MPICXX = CC
endif

# OpenMP
ifeq ($(COMPILERS), gcc-openmpi)
 OPENMP = -fopenmp
endif

ifeq ($(COMPILERS), gcc-sgi)
 OPENMP = -fopenmp
endif

ifeq ($(COMPILERS), gcc-cray)
 OPENMP = -fopenmp
endif

ifeq ($(COMPILERS), intel-intel)
 OPENMP = -qopenmp -mt_mpi
endif

ifeq ($(COMPILERS), intel-sgi)
 OPENMP = -qopenmp -mt_mpi
endif

ifeq ($(COMPILERS), intel-cray)
 OPENMP = -qopenmp
endif

# optimize or debug
OPTIMIZE = -O3 -DNDEBUG
#OPTIMIZE = -Wall -g
#gcc: full level of IEEE 754 compliance
#OPTIMIZE = -O3 -DNDEBUG -frounding-math -fsignaling-nans
#intel: -fp-model source also implies keyword precise; and it is nearly as twice slow as default -fp-model fast=1
#OPTIMIZE = -O3 -DNDEBUG -fp-model source

# GCC quadmath, default off
#QUADMATH = -DQUADMATH -L/usr/local/gcc-4.6.2/lib64 -lquadmath

# CXXFLAGS
CXXFLAGS = $(OPTIMIZE) $(OPENMP) $(QUADMATH)

###2. BOOST details
ifeq ($(PLATFORM), soilblast)
  ifeq ($(BOOST_ROOT), default)
   BOOST_LIB_DIR=-L/usr/lib64/openmpi/lib:/usr/lib64/compat-openmpi/lib:/usr/lib64
   BOOST_LIBS=-lboost_mpi -lboost_serialization -lboost_timer -lboost_system -lboost_chrono
   BOOST_RUN_LIB_DIR=-Wl,-rpath=/usr/lib64/openmpi/lib:/usr/lib64/compat-openmpi/lib:/usr/lib64
  else
   BOOST_INCLUDE=-I$(BOOST_ROOT)/include
   BOOST_LIB_DIR=-L$(BOOST_ROOT)/lib
   BOOST_LIBS=-lboost_mpi -lboost_serialization -lboost_timer -lboost_system -lboost_chrono
   BOOST_RUN_LIB_DIR=-Wl,-rpath=$(BOOST_ROOT)/lib

   ifeq ($(PAPI), yes)
    PAPI_INCLUDE=-I$(PAPI_ROOT)/include
    PAPI_LIB_DIR=-L$(PAPI_ROOT)/lib
    PAPI_LIBS=-lpapi
    PAPI_RUN_LIB_DIR=-Wl,-rpath=$(PAPI_ROOT)/lib
   endif
  endif
else
# for janus, spirit, excalibur, topaz, etc
 BOOST_INCLUDE=-I$(BOOST_ROOT)/include
 BOOST_LIB_DIR=-L$(BOOST_ROOT)/lib
 BOOST_LIBS=-lboost_mpi -lboost_serialization -lboost_timer -lboost_system -lboost_chrono
 BOOST_RUN_LIB_DIR=-Wl,-rpath=$(BOOST_ROOT)/lib

 ifeq ($(PAPI), yes)
  PAPI_INCLUDE=-I$(PAPI_ROOT)/include
  PAPI_LIB_DIR=-L$(PAPI_ROOT)/lib
  PAPI_LIBS=-lpapi
  PAPI_RUN_LIB_DIR=-Wl,-rpath=$(PAPI_ROOT)/lib
 endif
endif

INCLUDE=$(BOOST_INCLUDE) $(PAPI_INCLUDE)
LIB_DIR=$(BOOST_LIB_DIR) $(PAPI_LIB_DIR)
LIBS=$(BOOST_LIBS) $(PAPI_LIBS)
RUN_LIB_DIR=$(BOOST_RUN_LIB_DIR) $(PAPI_RUN_LIB_DIR)

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
