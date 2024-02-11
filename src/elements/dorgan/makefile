# $Id: makefile,v 1.13 2005-08-12 05:38:06 paklein Exp $
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
	J2SSKStV1D.o\
	J2_C0Hardening1DT.o\
	J2SSC0Hardening1DT.o\
	GradJ2SSKStV1D.o\
	GradJ2SSKStV2D.o\
	GradJ2SSKStV.o\
	GradSmallStrainT.o\
	GradSSSolidMatT.o\
	GradSSMatSupportT.o\
	GradSSSolidMatList1DT.o\
	GradSSSolidMatList2DT.o\
	GradSSSolidMatList3DT.o

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
	J2SSKStV1D.cpp\
	J2_C0Hardening1DT.cpp\
	J2SSC0Hardening1DT.cpp\
	GradJ2SSKStV1D.cpp\
	GradJ2SSKStV2D.cpp\
	GradJ2SSKStV.cpp\
	GradSmallStrainT.cpp\
	GradSSSolidMatT.cpp\
	GradSSMatSupportT.cpp\
	GradSSSolidMatList1DT.cpp\
	GradSSSolidMatList2DT.cpp\
	GradSSSolidMatList3DT.cpp

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
	J2SSKStV1D.h\
	J2_C0Hardening1DT.h\
	J2SSC0Hardening1DT.h\
	GradJ2SSKStV1D.h\
	GradJ2SSKStV2D.h\
	GradJ2SSKStV.h\
	GradSmallStrainT.h\
	GradSSSolidMatT.h\
	GradSSMatSupportT.h\
	GradSSSolidMatList1DT.h\
	GradSSSolidMatList2DT.h\
	GradSSSolidMatList3DT.h

HDR_LINK = $(HDR:.h=.h_link)

##################### (4) sub-directories
#
# (4.1) If there are sub-directories contained in this directory:
#    (1) list the subdirectories in the SUB_DIR macro
#    (2) uncomment the subdir_driver macro
#    (3) comment the SUB_DIR and subdir_driver in section (4.2)
#
SUB_DIR = \
        shape
subdir_driver: subdir_loop

#
# (4.2) If there are no sub-directories in this directory:
#    (1) comment the SUB_DIR and subdir_driver in section (4.1)
#    (2) uncomment the subdir_driver macro below
#
#SUB_DIR =
#subdir_driver: 

#####################

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = $(MACRO_DIR)/empty
include $(DEPEND)
