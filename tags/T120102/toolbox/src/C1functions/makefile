# $Id: makefile,v 1.22 2002-10-04 20:51:05 thao Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	CubicSplineT.o \
	C1FunctionT.o \
	ErrorFunc.o \
	GaoKlein.o \
	Gamma.o \
	GreenwoodWilliamson.o \
	ConHypGeom.o \
	LennardJones612.o \
	ModBessel.o \
	ParabolaT.o \
	PhiSplineT.o \
	SmithFerrante.o \
	ModSmithFerrante.o \
	LinearExponentialT.o \
	Triantafyllidis.o \
	LinearT.o \
	GaoVicky.o \
	GaoJi.o \
	GaoJi2.o \
	SF2.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	CubicSplineT.cpp \
	C1FunctionT.cpp \
	ErrorFunc.cpp \
	Gamma.cpp \
	GaoKlein.cpp \
	GreenwoodWilliamson.cpp \
	ConHypGeom.cpp \
	LennardJones612.cpp \
	ModBessel.cpp \
	ParabolaT.cpp \
	PhiSplineT.cpp \
	SmithFerrante.cpp \
	ModSmithFerrante.cpp \
	LinearExponentialT.cpp \
	Triantafyllidis.cpp \
	LinearT.cpp \
	GaoVicky.cpp \
	GaoJi.cpp \
	GaoJi2.cpp \
	SF2.cpp
# headers
HDR = \
	CubicSplineT.h \
	C1FunctionT.h \
	ErrorFunc.h \
	Gamma.h \
	GaoKlein.h \
	GreenwoodWilliamson.h \
	ConHypGeom.h \
	LennardJones612.h \
	ModBessel.h \
	ParabolaT.h \
	PhiSplineT.h \
	SmithFerrante.h \
	ModSmithFerrante.h \
	LinearExponentialT.h \
	Triantafyllidis.h \
	LinearT.h \
	GaoVicky.h \
	GaoJi.h \
	GaoJi2.h \
	SF2.h

HDR_LINK = $(HDR:.h=.h_link)

# subdirectories
SUB_DIR =  
#subdir_driver: subdir_loop
#SUB_DIR is empty
subdir_driver: 

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = /dev/null
include $(DEPEND)
