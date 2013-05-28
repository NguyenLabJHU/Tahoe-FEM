# options available in this makefile:
# 1. option to toggle on/off GCC quadmath (default off)
# 2. option to switch between optimization/debug mode (default optimization)

#PLATFORM = janus
PLATFORM = soilblast

ifeq ($(PLATFORM), janus)
 BOOST = boost-1.50.0_openmpi-1.4.3
endif

ifeq ($(PLATFORM), soilblast)
 BOOST = boost-1.53.0-openmpi-1.6.4-gcc-4.6.4
 #BOOST = boost-1.53.0
 #BOOST = default
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

 BOOST_ROOT=/home/yanb/local/$(BOOST)
 BOOST_INCLUDE=-I$(BOOST_ROOT)/include
 BOOST_LIB_DIR=-L$(BOOST_ROOT)/lib
 BOOST_LIBS=-lboost_mpi -lboost_serialization -lboost_timer
 BOOST_RUN_LIB_DIR=-Wl,-rpath=$(BOOST_ROOT)/lib

endif

ifeq ($(PLATFORM), soilblast)

  ifeq ($(BOOST), default)
   BOOST_LIB_DIR=-L/usr/lib64/openmpi/lib:/usr/lib64/compat-openmpi/lib:/usr/lib64
   BOOST_LIBS=-lboost_mpi -lboost_serialization -lboost_timer
   BOOST_RUN_LIB_DIR=-Wl,-rpath=/usr/lib64/openmpi/lib:/usr/lib64/compat-openmpi/lib:/usr/lib64
  else  
   BOOST_ROOT=/usr/local/$(BOOST)
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

# 4. general
SOURCES = $(wildcard *.cpp)
OBJECTS = $(SOURCES:.cpp=.o)
#OBJECTS = $(patsubst %.cpp, %.o, $(SOURCES))
EXECUTABLE = paraEllip3d

.PHONY: all tar clean

# template classes contain implementation and need to recompile upon change
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
