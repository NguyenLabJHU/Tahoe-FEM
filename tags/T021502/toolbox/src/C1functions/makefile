# $Id: makefile,v 1.11 2002-02-05 15:26:56 dzeigle Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	CubicSplineT.o \
	C1FunctionT.o \
	ConHyperGeom.o \
	ErrorFunc.o \
	GaoKlein.o \
	GenHyperGeom.o \
	GreenwoodWilliamson.o \
	LennardJones612.o \
	ModBessel.o \
	ParabolaT.o \
	PhiSplineT.o \
	SmithFerrante.o \
	LinearExponentialT.o \
	Triantafyllidis.o \
	LinearDecreaseT.o \
	QuadraticDecreaseT.o 

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	CubicSplineT.cpp \
	C1FunctionT.cpp \
	CONHyperGeom.cpp \
	ErrorFunc.cpp \
	GaoKlein.cpp \
	GenHyperGeom.cpp \
	GreenwoodWilliamson.cpp \
	LennardJones612.cpp \
	ModBessel.cpp \
	ParabolaT.cpp \
	PhiSplineT.cpp \
	SmithFerrante.cpp \
	LinearExponentialT.cpp \
	Triantafyllidis.cpp \
	LinearDecreaseT.cpp \
	QuadraticDecreaseT.cpp 

# headers
HDR = \
	CubicSplineT.h \
	C1FunctionT.h \
	ConHyperGeom.h \
	ErrorFunc.h \
	GaoKlein.h \
	GenHyperGeom.h \
	GreenwoodWilliamson.h \
	LennardJones612.h \
	ModBessel.h \
	ParabolaT.h \
	PhiSplineT.h \
	SmithFerrante.h \
	LinearExponentialT.h \
	Triantafyllidis.h \
	LinearDecreaseT.h \
	QuadraticDecreaseT.h 

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
