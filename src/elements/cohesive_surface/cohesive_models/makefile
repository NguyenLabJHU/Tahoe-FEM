# $Id: makefile,v 1.6 2001-10-11 23:19:52 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	LinearDamageT.o \
	SurfacePotentialT.o \
	TvergHutch2DT.o \
	ViscTvergHutch2DT.o \
	XuNeedleman2DT.o \
	XuNeedleman3DT.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	LinearDamageT.cpp \
	SurfacePotentialT.cpp \
	TvergHutch2DT.cpp \
	ViscTvergHutch2DT.cpp \
	XuNeedleman2DT.cpp \
	XuNeedleman3DT.cpp

# headers
HDR = \
	LinearDamageT.h \
	SurfacePotentialT.h \
	TvergHutch2DT.h \
	ViscTvergHutch2DT.h \
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
