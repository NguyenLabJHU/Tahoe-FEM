# $Id: makefile,v 1.8 2010-09-29 01:47:30 hspark Exp $
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
	CB_TersoffT.o \
	CB_TersoffT_surf.o \
	TersoffSolverT.o \
	TersoffSolverT_surf.o \
	dUdC.o \
	dUdC_surf.o \
	dXsi.o \
	dXsi_surf.o \
	ddC_driver.o \
	ddC_1.o \
	ddC_2.o \
	ddC_surf_driver.o \
	ddC_surf_1.o \
	ddC_surf_2.o \
	ddC_surf_3.o \
	ddC_surf_4.o \
	energy.o \
	surf_energy.o \
	pdm_energy.o \
	pdm_dUdC.o \
	pdm_dXsi.o \
	pdm_ddC.o

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
	CB_TersoffT.cpp \
	CB_TersoffT_surf.cpp \
	TersoffSolverT.cpp \
	TersoffSolverT_surf.cpp \
	dUdC.c \
	dUdC_surf.c \
	dXsi.c \
	dXsi_surf.c \
	ddC_driver.c \
	ddC_1.c \
	ddC_2.c \
	ddC_surf_driver.c \
	ddC_surf_1.c \
	ddC_surf_2.c \
	ddC_surf_3.c \
	ddC_surf_4.c \
	energy.c \
	surf_energy.c \
	pdm_energy.c \
	pdm_dUdC.c \
	pdm_dXsi.c \
	pdm_ddC.c

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
	CB_TersoffT.h \
	CB_TersoffT_surf.h \
	TersoffSolverT.h \
	TersoffSolverT_surf.h \
	Tersoff_common_defines_surf.h \
	Tersoff_common_defines.h \
	Tersoff_inc.h \
	Tersoff_inc_surf.h \
	PDM_common_defines.h \
	PDM_inc.h

HDR_LINK = $(HDR:.h=.h_link)

##################### (4) sub-directories
#
# (4.1) If there are sub-directories contained in this directory:
#    (1) list the subdirectories in the SUB_DIR macro
#    (2) uncomment the subdir_driver macro
#    (3) comment the SUB_DIR and subdir_driver in section (4.2)
#
SUB_DIR = \

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
