# options available in this makefile:
# 1. option to toggle on/off GCC quadmath (default off)
# 2. option to switch between optimization/debug mode (default optimization)

#PLATFORM = janus
PLATFORM = soilblast

ifeq ($(PLATFORM), janus)
 OPENMPI = 1.4.3
 #OPENMPI = 1.6
endif

# C++ compiler
MPICXX = mpic++

# optimization
OPTIMIZE = -O3 -DNDEBUG

# debugging options
DEBUG = -Wall -g

# OpenMP
OPENMP = -fopenmp

# 1. GCC quadmath
#QUADMATH = -DQUADMATH -L/usr/local/gcc-4.6.2/lib64 -lquadmath

# 2. Optimize or debug
CXXFLAGS = $(OPTIMIZE) $(OPENMP) $(QUADMATH)
#CXXFLAGS = $(DEBUG) $(OPENMP) $(QUADMATH)

# 3. Boost
ifeq ($(PLATFORM), janus)

 ifeq ($(OPENMPI), 1.4.3)

  BOOST_ROOT=/home/yanb/local/boost-1.49.0_openmpi-1.4.3
  BOOST_INCLUDE=-I$(BOOST_ROOT)/include
  BOOST_LIB_DIR=-L$(BOOST_ROOT)/lib -L/curc/tools/free/redhat_5_x86_64/openmpi-1.4.3_gcc-4.5.2_torque-2.5.8_ib/lib
  BOOST_LIBS=-lboost_mpi -lboost_serialization
  BOOST_RUN_LIB_DIR=-Wl,-rpath=$(BOOST_ROOT)/lib:/curc/tools/free/redhat_5_x86_64/openmpi-1.4.3_gcc-4.5.2_torque-2.5.8_ib/lib

 else 
  ifeq ($(OPENMPI), 1.6)

  BOOST_ROOT=/home/yanb/local/boost-1.49.0_openmpi-1.6
  BOOST_INCLUDE=-I$(BOOST_ROOT)/include
  BOOST_LIB_DIR=-L$(BOOST_ROOT)/lib -L/curc/tools/x_86_64/rh5/openmpi/1.6/gcc/4.7.1/torque/2.5.11/ib/lib
  BOOST_LIBS=-lboost_mpi -lboost_serialization
  BOOST_RUN_LIB_DIR=-Wl,-rpath=$(BOOST_ROOT)/lib:/curc/tools/x_86_64/rh5/openmpi/1.6/gcc/4.7.1/torque/2.5.11/ib/lib

  endif

 endif

else
 ifeq ($(PLATFORM), soilblast)

  BOOST_LIB_DIR=-L/usr/lib64/openmpi/lib:/usr/lib64/compat-openmpi/lib:/usr/lib64
  BOOST_LIBS=-lboost_mpi -lboost_serialization
  BOOST_RUN_LIB_DIR=-Wl,-rpath=/usr/lib64/openmpi/lib:/usr/lib64/compat-openmpi/lib:/usr/lib64
 endif

endif

INCLUDE=$(BOOST_INCLUDE)
LIB_DIR=$(BOOST_LIB_DIR)
LIBS=$(BOOST_LIBS)
RUN_LIB_DIR=$(BOOST_RUN_LIB_DIR)

# 4. general
SOURCES = $(wildcard *.cpp)
OBJECTS = $(SOURCES:.cpp=.o)
#OBJECTS = $(patsubst %.cpp, %.o, $(SOURCES))
EXECUTABLE = paraEllip3d

.PHONY: all tar clean

# template classes contain implementation and need to recompile upon change
all: $(EXECUTABLE) Contact.h Boundary.h

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
