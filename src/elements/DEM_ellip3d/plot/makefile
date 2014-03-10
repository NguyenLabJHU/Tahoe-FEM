BOOST_ROOT=/usr/local/boost-1.53.0_openmpi-1.6.4_gcc-4.6.4
#BOOST_ROOT=/home/yanb/local/boost-1.53.0_openmpi-1.6.4_gcc-4.7.2_ib

BOOST_INCLUDE=-I$(BOOST_ROOT)/include
BOOST_LIB_DIR=-L$(BOOST_ROOT)/lib
BOOST_RUN_LIB_DIR=-Wl,-rpath=$(BOOST_ROOT)/lib

INCLUDE=$(BOOST_INCLUDE)
LIB_DIR=$(BOOST_LIB_DIR)
RUN_LIB_DIR=$(BOOST_RUN_LIB_DIR)

MERGESRC = mergeContact.cpp
MERGEEXE = $(MERGESRC:.cpp=)

UNIQUESRC = uniqueContact.cpp
UNIQUEEXE = $(UNIQUESRC:.cpp=)

TOOLSRC = $(filter-out $(MERGESRC) $(UNIQUESRC), $(wildcard *.cpp) )
TOOLEXE = $(TOOLSRC:.cpp=)

CXX = g++
CXXFLAGS =

.PHONY: all clean

all: $(TOOLEXE)

%: %.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

mergeContact: $(MERGEEXE)
uniqueContact: $(UNIQUEEXE)

$(MERGEEXE): $(MERGESRC)
	$(CXX) $(CXXFLAGS) $< -o $@ $(INCLUDE) $(LIB_DIR) $(RUN_LIB_DIR)

$(UNIQUEEXE): $(UNIQUESRC)
	$(CXX) $(CXXFLAGS) $< -o $@ $(INCLUDE) $(LIB_DIR) $(RUN_LIB_DIR)

clean:
	-rm -f *~ $(TOOLEXE) $(MERGEEXE) $(UNIQUEEXE)
