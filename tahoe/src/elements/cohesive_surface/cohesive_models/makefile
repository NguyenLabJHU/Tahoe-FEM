# $Id: makefile,v 1.18 2004-07-15 08:26:03 paklein Exp $
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
	RateDep2DT.o \
	LinearDamageT.o \
	SurfacePotentialT.o \
	SurfacePotentialT_factory.o \
	TiedPotentialT.o \
	Tijssens2DT.o \
	TvergHutch2DT.o \
	ViscTvergHutch2DT.o \
	XuNeedleman2DT.o \
	XuNeedleman3DT.o \
	YoonAllen2DT.o \
	YoonAllen3DT.o \
	TvergHutch3DT.o \
	TiedPotentialBaseT.o \
	From2Dto3DT.o \
	TvergHutchRigid2DT.o

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)

#
# (1.2) If there are not object in this directory:
#    (1) comment the OBJ, OBJ_LINK, and DEP definitions in section (1.1)
#    (2) uncomment the OBJ, OBJ_LINK, and DEP definitions below
#
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

##################### (2) sources
#
# (2) If there are sources in this directory list them in the
#     SRC macro; otherwise, comment the SRC macro.
#
SRC = \
	RateDep2DT.cpp \
	LinearDamageT.cpp \
	SurfacePotentialT.cpp \
	SurfacePotentialT_factory.cpp \
	TiedPotentialT.cpp \
	Tijssens2DT.cpp \
	TvergHutch2DT.cpp \
	ViscTvergHutch2DT.cpp \
	XuNeedleman2DT.cpp \
	XuNeedleman3DT.cpp \
	YoonAllen2DT.cpp \
	YoonAllen3DT.cpp \
	TvergHutch3DT.cpp \
	TiedPotentialBaseT.cpp \
	From2Dto3DT.cpp \
	TvergHutchRigid2DT.cpp

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
	TvergHutch3DT.h \
	TiedPotentialBaseT.h \
	From2Dto3DT.h \
	TvergHutchRigid2DT.h

HDR_LINK = $(HDR:.h=.h_link)

##################### (4) sub-directories
#
# (4.1) If there are sub-directories contained in this directory:
#    (1) list the subdirectories in the SUB_DIR macro
#    (2) uncomment the subdir_driver macro
#    (3) comment the SUB_DIR and subdir_driver in section (4.2)
#
#SUB_DIR = 
#subdir_driver: subdir_loop

#
# (4.2) If there are no sub-directories in this directory:
#    (1) comment the SUB_DIR and subdir_driver in section (4.1)
#    (2) uncomment the subdir_driver macro below
#
SUB_DIR =
subdir_driver: 

#####################

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = /dev/null
include $(DEPEND)
