# $Id: makefile,v 1.4 2001-02-22 23:46:17 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	DecohesionLawT.o \
	LinearDamageT.o \
	SurfacePotentialT.o \
	TvergHutch2DT.o \
	XuNeedleman2DT.o \
	XuNeedleman3DT.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	DecohesionLawT.cpp \
	LinearDamageT.cpp \
	SurfacePotentialT.cpp \
	TvergHutch2DT.cpp \
	XuNeedleman2DT.cpp \
	XuNeedleman3DT.cpp

# headers
HDR = \
	DecohesionLawT.h \
	LinearDamageT.h \
	SurfacePotentialT.h \
	TvergHutch2DT.h \
	XuNeedleman2DT.h \
	XuNeedleman3DT.h
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
