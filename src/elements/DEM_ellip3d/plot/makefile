#Usage: 
#  make
#  make merge
#  make unique

###1. system platform
PLATFORM = soilblast
#PLATFORM = spirit
#PLATFORM = janus

###2. compilers package
COMPILERS = gcc
#COMPILERS = intel

ifeq ($(COMPILERS), gcc)
 CXX = g++
endif

ifeq ($(COMPILERS), intel)
 CXX = icpc
endif

CXXFLAGS =

###3. boost location
ifeq ($(PLATFORM), soilblast)
 BOOST_ROOT = /usr/local/boost-1.53.0-openmpi-1.6.4-gcc-4.6.4
 #BOOST_ROOT = default
endif

ifeq ($(PLATFORM), spirit)
 ifeq ($(COMPILERS), gcc)
  BOOST_ROOT = /work4/projects/openmpi/boost-1.55.0-openmpi-1.7.5-gcc-compilers-4.7.3
 endif
 ifeq ($(COMPILERS), intel)
  #BOOST_ROOT = /work4/projects/openmpi/boost-1.53.0-intelmpi-4.0.3.008
  BOOST_ROOT = /work4/projects/openmpi/boost-1.55.0-intelmpi-4.1.3.048-intel-compilers-13.0.1
 endif
endif

ifeq ($(PLATFORM), janus)
 ifeq ($(COMPILERS), gcc)
  #BOOST_ROOT = /home/yanb/local/boost-1.55.0_openmpi-1.7.3_gcc-4.8.2_ib
  BOOST_ROOT = /home/yanb/local/boost-1.55.0_openmpi-1.8.2_gcc-4.9.1
 endif
 ifeq ($(COMPILERS), intel)
  BOOST_ROOT = /home/yanb/local/boost-1.55.0_impi-13.0.0_intel-13.0.0
 endif
endif

BOOST_INCLUDE=-I$(BOOST_ROOT)/include
BOOST_LIB_DIR=-L$(BOOST_ROOT)/lib
BOOST_RUN_LIB_DIR=-Wl,-rpath=$(BOOST_ROOT)/lib

###4. makefile
INCLUDE=$(BOOST_INCLUDE)
LIB_DIR=$(BOOST_LIB_DIR)
RUN_LIB_DIR=$(BOOST_RUN_LIB_DIR)

MERGESRC = mergeContact.cpp
MERGEEXE = $(MERGESRC:.cpp=)

UNIQUESRC = uniqueContact.cpp
UNIQUEEXE = $(UNIQUESRC:.cpp=)

TOOLSRC = $(filter-out $(MERGESRC) $(UNIQUESRC), $(wildcard *.cpp) )
TOOLEXE = $(TOOLSRC:.cpp=)

.PHONY: all clean

all: $(TOOLEXE)

%: %.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

merge: $(MERGEEXE)
unique: $(UNIQUEEXE)

$(MERGEEXE): $(MERGESRC)
	$(CXX) $(CXXFLAGS) $< -o $@ $(INCLUDE) $(LIB_DIR) $(RUN_LIB_DIR)

$(UNIQUEEXE): $(UNIQUESRC)
	$(CXX) $(CXXFLAGS) $< -o $@ $(INCLUDE) $(LIB_DIR) $(RUN_LIB_DIR)

clean:
	-rm -f *~ $(TOOLEXE) $(MERGEEXE) $(UNIQUEEXE)
