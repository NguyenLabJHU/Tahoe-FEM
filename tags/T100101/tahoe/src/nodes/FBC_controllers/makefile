# $Id: makefile,v 1.5 2001-09-11 23:43:48 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	AugLagSphereT.o \
	AugLagWallT.o \
	FBC_ControllerT.o \
	MFPenaltySphereT.o \
	PenaltyRegionT.o \
	PenaltySphereT.o \
	PenaltyWallT.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	AugLagSphereT.cpp \
	AugLagWallT.cpp \
	FBC_ControllerT.cpp \
	MFPenaltySphereT.cpp \
	PenaltyRegionT.cpp \
	PenaltySphereT.cpp \
	PenaltyWallT.cpp

# headers
HDR = \
	AugLagSphereT.h \
	AugLagWallT.h \
	FBC_ControllerT.h \
	MFPenaltySphereT.h \
	PenaltyRegionT.h \
	PenaltySphereT.h \
	PenaltyWallT.h
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
