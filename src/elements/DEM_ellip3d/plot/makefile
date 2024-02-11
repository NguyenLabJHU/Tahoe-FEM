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
 BOOST_ROOT = /usr/local/boost-1.70.0-openmpi-4.0.2-gcc-9.2.0
 #BOOST_ROOT = default
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
