#########################################################
###        Part A: change with your environment       ###
#########################################################
###1. system platform
PLATFORM = soilblast
#PLATFORM = topaz
#PLATFORM = onyx
#PLATFORM = excalibur
#PLATFORM = thunder
#PLATFORM = centennial
#PLATFORM = lightning
#PLATFORM = garnet

###2. compiler combinations (c++ and mpi)
COMPILERS = gcc-openmpi
#COMPILERS = intel-cray
#COMPILERS = intel-sgi
#COMPILERS = intel-intel
#COMPILERS = gcc-sgi
#COMPILERS = gcc-cray

###3. compute continuum_stress-strain or not?
STRESS_STRAIN = yes
ifeq ($(STRESS_STRAIN), yes)
 EIGEN = yes
 QHULL = yes
 STRESS_STRAIN_MACRO = -DSTRESS_STRAIN
endif

###4. use papi or not?
#PAPI = yes

###5. boost location
ifeq ($(PLATFORM), soilblast)
 ifeq ($(COMPILERS), gcc-openmpi)
  #BOOST_ROOT = default
  #BOOST_ROOT = /usr/local/boost-1.53.0-openmpi-1.6.4-gcc-4.6.4
  BOOST_ROOT = /usr/local/boost-1.65.1-openmpi-3.0.0-gcc-6.4.0
  ifeq ($(PAPI), yes)
   PAPI_ROOT = /usr/local/papi-5.5.1-openmpi-3.0.0-gcc-6.4.0
  endif
  ifeq ($(EIGEN), yes)
   EIGEN_ROOT = /usr/local/eigen-3.3.4
  endif
  ifeq ($(QHULL), yes)
   QHULL_ROOT = /usr/local/qhull-2015.2-gcc-6.4.0
  endif
 endif
endif

ifeq ($(PLATFORM), excalibur)
 ifeq ($(COMPILERS), intel-cray)
  BOOST_ROOT = /p/home/yanb/local/boost-1.65.1_cray-mpich7.1.0_PrgEnv-intel5.2.40
  ifeq ($(PAPI), yes)
   PAPI_ROOT = /p/home/yanb/local/papi-5.4.1_cray-mpich7.1.0_PrgEnv-intel5.2.40
  endif
  ifeq ($(EIGEN), yes)
   EIGEN_ROOT = /p/home/yanb/local/eigen-3.3.4
  endif
  ifeq ($(QHULL), yes)
   QHULL_ROOT = /p/home/yanb/local/qhull-2015.2-intel5.2.40
  endif
 endif
 ifeq ($(COMPILERS), gcc-cray)
  BOOST_ROOT = /p/home/yanb/local/boost-1.65.1_cray-mpich7.1.0_PrgEnv-gnu5.2.40
  ifeq ($(EIGEN), yes)
   EIGEN_ROOT = /p/home/yanb/local/eigen-3.3.4
  endif
  ifeq ($(QHULL), yes)
   QHULL_ROOT = /p/home/yanb/local/qhull-2015.2-gnu5.2.40
  endif
 endif
 ifeq ($(COMPILERS), gcc-openmpi)
  BOOST_ROOT = /p/home/yanb/local/boost-1.57.0_openmpi-1.8.4_PrgEnv-gnu5.2.40
  ifeq ($(PAPI), yes)
   PAPI_ROOT = /p/home/yanb/local/papi-5.5.1_openmpi-1.8.4_PrgEnv-gnu5.2.40 
  endif
  ifeq ($(EIGEN), yes)
   EIGEN_ROOT = /p/home/yanb/local/eigen-3.3.4
  endif
  ifeq ($(QHULL), yes)
   QHULL_ROOT = /p/home/yanb/local/qhull-2015.2-gnu5.2.40
  endif
 endif
endif

ifeq ($(PLATFORM), onyx)
 ifeq ($(COMPILERS), intel-cray)
  BOOST_ROOT = /p/home/yanb/local/boost-1.65.1_cray-mpich7.6.2_PrgEnv-intel6.0.4
  ifeq ($(EIGEN), yes)
   EIGEN_ROOT = /p/home/yanb/local/eigen-3.3.4
  endif
  ifeq ($(QHULL), yes)
   QHULL_ROOT = /p/home/yanb/local/qhull-2015.2-intel6.0.4
  endif
 endif
 ifeq ($(COMPILERS), gcc-cray)
  BOOST_ROOT = /p/home/yanb/local/boost-1.65.1_cray-mpich7.6.2_PrgEnv-gnu6.0.4
  ifeq ($(EIGEN), yes)
   EIGEN_ROOT = /p/home/yanb/local/eigen-3.3.4
  endif
  ifeq ($(QHULL), yes)
   QHULL_ROOT = /p/home/yanb/local/qhull-2015.2-gnu6.0.4
  endif
 endif
 ifeq ($(COMPILERS), gcc-openmpi)
  BOOST_ROOT = 
  ifeq ($(PAPI), yes)
   PAPI_ROOT = 
  endif
 endif
endif

ifeq ($(PLATFORM), topaz)
 ifeq ($(COMPILERS), intel-sgi)
  BOOST_ROOT = /p/home/yanb/local/boost-1.65.1_mpt-2.14_intel-16.0.0
  ifeq ($(PAPI), yes)
   PAPI_ROOT = /p/home/yanb/local/papi-5.4.1_mpt-2.14_intel-16.0.0
  endif
  ifeq ($(EIGEN), yes)
   EIGEN_ROOT = /p/home/yanb/local/eigen-3.3.4
  endif
  ifeq ($(QHULL), yes)
   QHULL_ROOT = /p/home/yanb/local/qhull-2015.2-intel-16.0.0
  endif
 endif
endif

ifeq ($(PLATFORM), thunder)
 ifeq ($(COMPILERS), intel-sgi)
  BOOST_ROOT = /home/yanb/local/boost-1.65.1_mpt-2.14_intel-compilers-15.0.3.187
  ifeq ($(PAPI), yes)
   PAPI_ROOT = /home/yanb/local/papi-5.4.1_mpt-2.14_intel-compilers-15.0.3.187
  endif
  ifeq ($(EIGEN), yes)
   EIGEN_ROOT = /home/yanb/local/eigen-3.3.4
  endif
  ifeq ($(QHULL), yes)
   QHULL_ROOT = /home/yanb/local/qhull-2015.2-intel-compilers-15.0.3.187
  endif
 endif
 ifeq ($(COMPILERS), gcc-sgi)
  BOOST_ROOT = /home/yanb/local/boost-1.57.0_mpt-2.12_gcc-compilers-4.4.6
 endif
endif

ifeq ($(PLATFORM), centennial)
 ifeq ($(COMPILERS), intel-sgi)
  BOOST_ROOT = /p/home/yanb/local/boost-1.65.1_mpt-2.15_intel-2017.1.132
  ifeq ($(EIGEN), yes)
   EIGEN_ROOT = /p/home/yanb/local/eigen-3.3.4
  endif
  ifeq ($(QHULL), yes)
   QHULL_ROOT = /p/home/yanb/local/qhull-2015.2-intel-2017.1.132
  endif
 endif
endif

ifeq ($(PLATFORM), lightning)
 ifeq ($(COMPILERS), intel-cray)
  #BOOST_ROOT = /home/yanb/local/boost-1.55.0_cray-mpich7.1.3_PrgEnv-intel5.2.40
  BOOST_ROOT = /home/yanb/local/boost-1.57.0_cray-mpich7.2.6_PrgEnv-intel5.2.82
 endif
endif

ifeq ($(PLATFORM), garnet)
 ifeq ($(COMPILERS), intel-cray)
  BOOST_ROOT = /u/yanb/local/boost-1.57.0_cray-mpich7.1.0_PrgEnv-intel5.2.82
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

###2 standard, optimize or debug
#C++11 or not?
C++11 = -std=c++11

OPTIMIZE = -O3 -DNDEBUG $(C++11)
#OPTIMIZE = -Wall -g $(C++11)
#gcc: full level of IEEE 754 compliance
#OPTIMIZE = -O3 -DNDEBUG -frounding-math -fsignaling-nans $(C++11)
#intel: -fp-model source also implies keyword precise; and it is nearly as twice slow as default -fp-model fast=1
#OPTIMIZE = -O3 -DNDEBUG -fp-model source $(C++11)

# GCC quadmath, default off
#QUADMATH = -DQUADMATH -L/usr/local/gcc-4.6.2/lib64 -lquadmath

# CXXFLAGS
CXXFLAGS = $(OPTIMIZE) $(OPENMP) $(QUADMATH) $(STRESS_STRAIN_MACRO)

###3. configurations of BOOST, PAPI, Eigen and Qhull
ifeq ($(PLATFORM), soilblast)
  ifeq ($(BOOST_ROOT), default)
   BOOST_LIB_DIR=-L/usr/lib64/openmpi/lib:/usr/lib64/compat-openmpi/lib:/usr/lib64
   BOOST_LIBS=-lboost_mpi -lboost_serialization -lboost_timer -lboost_chrono -lboost_system
   BOOST_RUN_LIB_DIR=-Wl,-rpath=/usr/lib64/openmpi/lib:/usr/lib64/compat-openmpi/lib:/usr/lib64
  else
   BOOST_INCLUDE=-I$(BOOST_ROOT)/include
   BOOST_LIB_DIR=-L$(BOOST_ROOT)/lib
   BOOST_LIBS=-lboost_mpi -lboost_serialization -lboost_timer -lboost_chrono -lboost_system
   BOOST_RUN_LIB_DIR=-Wl,-rpath=$(BOOST_ROOT)/lib

   ifeq ($(PAPI), yes)
    PAPI_INCLUDE=-I$(PAPI_ROOT)/include
    PAPI_LIB_DIR=-L$(PAPI_ROOT)/lib
    PAPI_LIBS=-lpapi
    PAPI_RUN_LIB_DIR=-Wl,-rpath=$(PAPI_ROOT)/lib
   endif

   ifeq ($(EIGEN), yes)
    EIGEN_INCLUDE=-I$(EIGEN_ROOT)
   endif

   ifeq ($(QHULL), yes)
    QHULL_INCLUDE=-I$(QHULL_ROOT)/src
    QHULL_LIB_DIR=-L$(QHULL_ROOT)/lib
    #QHULL_LIBS=-lqhullcpp -lqhullstatic_r
    QHULL_LIBS=-lqhullcpp -lqhull_r
    QHULL_RUN_LIB_DIR=-Wl,-rpath=$(QHULL_ROOT)/lib	
   endif

  endif
else
# for excalibur, topaz, etc
 BOOST_INCLUDE=-I$(BOOST_ROOT)/include
 BOOST_LIB_DIR=-L$(BOOST_ROOT)/lib
 BOOST_LIBS=-lboost_mpi -lboost_serialization -lboost_timer -lboost_chrono -lboost_system
 BOOST_RUN_LIB_DIR=-Wl,-rpath=$(BOOST_ROOT)/lib

 ifeq ($(PAPI), yes)
  PAPI_INCLUDE=-I$(PAPI_ROOT)/include
  PAPI_LIB_DIR=-L$(PAPI_ROOT)/lib
  PAPI_LIBS=-lpapi
  PAPI_RUN_LIB_DIR=-Wl,-rpath=$(PAPI_ROOT)/lib
 endif

 ifeq ($(EIGEN), yes)
  EIGEN_INCLUDE=-I$(EIGEN_ROOT)
 endif

 ifeq ($(QHULL), yes)
  QHULL_INCLUDE=-I$(QHULL_ROOT)/src
  QHULL_LIB_DIR=-L$(QHULL_ROOT)/lib
  ## it depends on the type of Boost libraries; or use -Bstatic -Bdynamic to control.
  QHULL_LIBS=-lqhullcpp -lqhullstatic_r
  #QHULL_LIBS=-lqhullcpp -lqhull_r
  #QHULL_RUN_LIB_DIR=-Wl,-rpath=$(QHULL_ROOT)/lib
 endif

endif

INCLUDE=$(BOOST_INCLUDE) $(PAPI_INCLUDE) $(EIGEN_INCLUDE) $(QHULL_INCLUDE)
LIB_DIR=$(BOOST_LIB_DIR) $(PAPI_LIB_DIR) $(QHULL_LIB_DIR)
LIBS=$(BOOST_LIBS) $(PAPI_LIBS) $(QHULL_LIBS)
RUN_LIB_DIR=$(BOOST_RUN_LIB_DIR) $(PAPI_RUN_LIB_DIR) $(QHULL_RUN_LIB_DIR)

###4. makefile
SOURCES = $(wildcard *.cpp)
OBJECTS = $(SOURCES:.cpp=.o)
#OBJECTS = $(patsubst %.cpp, %.o, $(SOURCES))
EXECUTABLE = paraEllip3d

.PHONY: all tar clean

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(MPICXX) -o $@ $(CXXFLAGS) $(OBJECTS) $(LIB_DIR) $(LIBS) $(RUN_LIB_DIR)

%.o: %.cpp
	$(MPICXX) -c $< -o $@ $(CXXFLAGS) $(INCLUDE)

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
