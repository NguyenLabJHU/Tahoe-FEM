# $Id: makefile,v 1.6 2003-01-28 22:35:26 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	ThermalSurfaceT.o \
	CSEAnisoT.o \
	CSEBaseT.o \
	CSEIsoT.o \
	MeshFreeCSEAnisoT.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	ThermalSurfaceT.cpp \
	CSEAnisoT.cpp \
	CSEBaseT.cpp \
	CSEIsoT.cpp \
	MeshFreeCSEAnisoT.cpp

# headers
HDR = \
	ThermalSurfaceT.h \
	CSEAnisoT.h \
	CSEBaseT.h \
	CSEIsoT.h \
	MeshFreeCSEAnisoT.h
HDR_LINK = $(HDR:.h=.h_link)

# subdirectories
SUB_DIR = \
	cohesive_models

subdir_driver: subdir_loop
# SUB_DIR is empty
#subdir_driver: 

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = /dev/null
include $(DEPEND)
