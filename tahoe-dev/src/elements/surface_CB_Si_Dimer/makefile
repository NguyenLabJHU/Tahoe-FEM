# $Id: makefile,v 1.4 2008-09-02 13:34:33 hspark Exp $
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
	CB_TersoffDimerT.o \
	CB_TersoffDimerT_surf.o \
	TersoffDimerSolverT.o \
	TersoffDimerSolverT_surf.o \
	TDdUdC.o \
	TDdUdC_surf.o \
	TDdXsi.o \
	TDdXsi_surf_driver.o \
	TDdXsi_surf_1.o \
	TDdXsi_surf_2.o \
	TDdXsi_surf_3.o \
	TDdXsi_surf_4.o \
	TDdXsi_surf_5.o \
	TDdXsi_surf_6.o \
	TDdXsi_surf_7.o \
	TDdXsi_surf_8.o \
	TDddC.o \
	TDddC_surf_driver.o \
	TDddC_surf_1.o \
	TDddC_surf_2.o \
	TDddC_surf_3.o \
	TDddC_surf_4.o \
	TDddC_surf_5.o \
	TDddC_surf_6.o \
	TDddC_surf_7.o \
	TDddC_surf_8.o \
	TDddC_surf_9.o \
	TDddC_surf_10.o \
	TDddC_surf_11.o \
	TDddC_surf_12.o \
	TDddC_surf_13.o \
	TDddC_surf_14.o \
	TDddC_surf_15.o \
	si_dimer_energy.o 

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
	CB_TersoffDimerT.cpp \
	CB_TersoffDimerT_surf.cpp \
	TersoffSolverDimerT.cpp \
	TersoffSolverDimerT_surf.cpp \
	TDdUdC.c \
	TDdUdC_surf.c \
	TDdXsi.c \
	TDdXsi_surf_driver.c \
	TDdXsi_surf_1.c \
	TDdXsi_surf_2.c \
	TDdXsi_surf_3.c \
	TDdXsi_surf_4.c \
	TDdXsi_surf_5.c \
	TDdXsi_surf_6.c \
	TDdXsi_surf_7.c \
	TDdXsi_surf_8.c \
	TDddC.c \
	TDddC_surf_driver.c \
	TDddC_surf_1.c \
	TDddC_surf_2.c \
	TDddC_surf_3.c \
	TDddC_surf_4.c \
	TDddC_surf_5.c \
	TDddC_surf_6.c \
	TDddC_surf_7.c \
	TDddC_surf_8.c \
	TDddC_surf_9.c \
	TDddC_surf_10.c \
	TDddC_surf_11.c \
	TDddC_surf_12.c \
	TDddC_surf_13.c \
	TDddC_surf_14.c \
	TDddC_surf_15.c \
	si_dimer_energy.c


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
	CB_TersoffDimerT.h \
	CB_TersoffDimerT_surf.h \
	TersoffDimerSolverT.h \
	TersoffDimerSolverT_surf.h \
	TersoffDimer_common_defines_surf.h \
	TersoffDimer_common_defines.h \
	TersoffDimer_inc.h \
	TersoffDimer_inc_surf.h

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
