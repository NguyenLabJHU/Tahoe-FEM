# $Id: makefile,v 1.32 2009-05-20 17:48:26 regueiro Exp $
#
# This makefile is divided into four sections:
#    (1) objects
#    (2) sources
#    (3) headers
#    (4) subdirectories
#
# Instructions for how to define these is given in each section below
#
include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

##################### (1) objects
#
# (1.1) If there are objects in this directory:
#    (1) define the OBJ macro as the list of object files
#    (2) uncomment the OBJ_LINK and DEP macros here
#    (3) comment the OBJ, OBJ_LINK, and DEP macros in the
#        section (1.1).
#
OBJ = \
	CubicSplineT.o \
	C1FunctionT.o \
	C1FunctionT_factory.o \
	ErrorFunc.o \
	GaoKlein.o \
	Gamma.o \
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
	SF2.o \
	PiecewiseLinearT.o \
	PowerLawT.o\
	CosineT.o \
	CosinePlusT.o \
        InvLangevin.o

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)

#
# (1.2) If there are not object in this directory:
#    (1) comment the OBJ, OBJ_LINK, and DEP definitions in section (1.1)
#    (2) uncomment the OBJ, OBJ_LINK, and DEP definitions below
#
#OBJ = dummy
#OBJ_LINK =
#DEP = $(MACRO_DIR)/empty

##################### (2) sources
#
# (2) If there are sources in this directory list them in the
#     SRC macro; otherwise, comment the SRC macro.
#
SRC = \
	CubicSplineT.cpp \
	C1FunctionT.cpp \
	C1FunctionT_factory.cpp \
	ErrorFunc.cpp \
	Gamma.cpp \
	GaoKlein.cpp \
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
	SF2.cpp \
	PiecewiseLinearT.cpp \
	PowerLawT.cpp \
	CosineT.cpp \
	CosinePlusT.cpp \
        InvLangevin.cpp

##################### (3) headers
#
# (3) If there are header files in this directory:
#    (1) list the header files in the HDR macro
#    (2) uncomment the HDR_LINK
# 
# If there are no header files in this directory comment both the
# HDR and HDR_LINK macro   
#
HDR = \
	CubicSplineT.h \
	C1FunctionT.h \
	ErrorFunc.h \
	Gamma.h \
	GaoKlein.h \
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
	SF2.h \
	PiecewiseLinearT.h \
	PowerLawT.h \
	CosineT.h \
	CosinePlusT.h \
        InvLangevin.h

HDR_LINK = $(HDR:.h=.h_link)

##################### (4) sub-directories
#
# (4.1) If there are sub-directories contained in this directory:
#    (1) list the subdirectories in the SUB_DIR macro
#    (2) uncomment the subdir_driver macro
#    (3) comment the SUB_DIR and subdir_driver in section (4.2)
#
SUB_DIR =  
#subdir_driver: subdir_loop

#
# (4.2) If there are no sub-directories in this directory:
#    (1) comment the SUB_DIR and subdir_driver in section (4.1)
#    (2) uncomment the subdir_driver macro below
#
#SUB_DIR is empty
subdir_driver: 

#####################

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = $(MACRO_DIR)/empty
include $(DEPEND)
