# $Id: makefile,v 1.7 2001-06-11 01:55:49 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	CubicSplineT.o \
	C1FunctionT.o \
	GaoKlein.o \
	LennardJones612.o \
	ParabolaT.o \
	PhiSplineT.o \
	SmithFerrante.o \
	LinearExponentialT.o \
	Triantafyllidis.o

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	CubicSplineT.cpp \
	C1FunctionT.cpp \
	GaoKlein.cpp \
	LennardJones612.cpp \
	ParabolaT.cpp \
	PhiSplineT.cpp \
	SmithFerrante.cpp \
	LinearExponentialT.cpp \
	Triantafyllidis.cpp

# headers
HDR = \
	CubicSplineT.h \
	C1FunctionT.h \
	GaoKlein.h \
	LennardJones612.h \
	ParabolaT.h \
	PhiSplineT.h \
	SmithFerrante.h \
	LinearExponentialT.h \
	Triantafyllidis.h

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
