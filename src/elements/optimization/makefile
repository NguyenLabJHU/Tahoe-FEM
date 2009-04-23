# $Header: /home/regueiro/tahoe_cloudforge_repo_snapshots/development/src/elements/optimization/makefile,v 1.1 2009-04-23 03:03:44 thao Exp $
# created: tdnguye (07/04/2006)
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
OBJ =  SS_Optimize_Dual.o \
	FSFiber_Optimize_Dual.o \
	FSFiber_OptSurf.o \
	FSFiber_OptNS.o 

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)

#OBJ = dummy
#OBJ_LINK =
#DEP = $(MACRO_DIR)/empty

##################### (2) sources
#
# (2) If there are sources in this directory list them in the
#     SRC macro; otherwise, comment the SRC macro.
#
SRC = 	SS_Optimize_Dual.cpp \
		FSFiber_Optimize_Dual.cpp \
		FSFiber_OptSurf.cpp \
		FSFiber_OptNS.cpp

##################### (3) headers
#
# (3) If there are header files in this directory:
#    (1) list the header files in the HDR macro
#    (2) uncomment the HDR_LINK
# 
# If there are no header files in this directory comment both the
# HDR and HDR_LINK macro   
#
HDR = 	SS_Optimize_Dual.h \
	FSFiber_Optimize_Dual.h \
	FSFiber_OptSurf.h \
	FSFiber_OptNS.h

HDR_LINK = $(HDR:.h=.h_link)

##################### (4) sub-directories
#
# (4.1) If there are sub-directories contained in this directory:
#    (1) list the subdirectories in the SUB_DIR macro
#    (2) uncomment the subdir_driver macro
#    (3) comment the SUB_DIR and subdir_driver in section (4.2)
#
SUB_DIR = \
	materials

subdir_driver: subdir_loop
#
# (4.2) If there are no sub-directories in this directory:
#    (1) comment the SUB_DIR and subdir_driver in section (4.1)
#    (2) uncomment the subdir_driver macro below
#
#SUB_DIR =
#subdir_driver: 
# SUB_DIR is empty
#subdir_driver: 

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = $(MACRO_DIR)/empty
include $(DEPEND)
