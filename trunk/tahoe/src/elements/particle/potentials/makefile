# $Id: makefile,v 1.1 2002-11-26 01:55:37 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	HarmonicPairT.o \
	LennardJonesPairT.o \
	PairPropertyT.o \
	ParticlePropertyT.o

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	HarmonicPairT.cpp \
	LennardJonesPairT.cpp \
	PairPropertyT.cpp \
	ParticlePropertyT.cpp

# headers
HDR = \
	HarmonicPairT.h \
	LennardJonesPairT.h \
	PairPropertyT.h \
	ParticlePropertyT.h

HDR_LINK = $(HDR:.h=.h_link)

# subdirectories
SUB_DIR = 
#subdir_driver: subdir_loop
# SUB_DIR is empty
subdir_driver: 

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = /dev/null
include $(DEPEND)
