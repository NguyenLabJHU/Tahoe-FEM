# $Id: makefile,v 1.13 2002-08-05 19:44:05 cjkimme Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	RateDep2DT.o \
	LinearDamageT.o \
	SurfacePotentialT.o \
	TiedPotentialT.o \
	Tijssens2DT.o \
	TvergHutch2DT.o \
	ViscTvergHutch2DT.o \
	XuNeedleman2DT.o \
	XuNeedleman3DT.o \
	YoonAllen2DT.o \
	YoonAllen3DT.o \
	TvergHutch3DT.o 

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	RateDep2DT.cpp \
	LinearDamageT.cpp \
	SurfacePotentialT.cpp \
	TiedPotentialT.cpp \
	Tijssens2DT.cpp \
	TvergHutch2DT.cpp \
	ViscTvergHutch2DT.cpp \
	XuNeedleman2DT.cpp \
	XuNeedleman3DT.cpp \
	YoonAllen2DT.cpp \
	YoonAllen3DT.cpp \
	TvergHutch3DT.cpp

# headers
HDR = \
	RateDep2DT.h \
	LinearDamageT.h \
	SurfacePotentialT.h \
	TiedPotentialT.h \
	Tijssens2DT.h \
	TvergHutch2DT.h \
	ViscTvergHutch2DT.h \
	XuNeedleman2DT.h \
	XuNeedleman3DT.h \
	YoonAllen2DT.h \
	YoonAllen3DT.h \
	TvergHutch3DT.h

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
