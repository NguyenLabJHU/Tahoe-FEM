# $Id: makefile,v 1.10 2004-07-22 08:32:55 paklein Exp $
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
	ExecutionManagerT.o \
	FEExecutionManagerT.o \
	FEManagerT.o \
	FEManagerT.ParseInput.o \
	FEManagerT_mpi.o \
	FEManagerT_bridging.o \
	FEManagerT_bridging_1.o \
	FEManagerT_bridging_2.o \
	FEManagerT_bridging_22.o \
	FEManagerT_bridging_3.o \
	FEManagerT_bridging_4.o \
	GlobalT.o \
	IOManager.o \
	IOManager_mpi.o \
	main.o \
	TimeManagerT.o \
	CommManagerT.o \
	MultiManagerT.o

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
	ExecutionManagerT.cpp \
	FEExecutionManagerT.cpp \
	FEManagerT.cpp \
	FEManagerT.ParseInput.cpp \
	FEManagerT_mpi.cpp \
	FEManagerT_bridging.cpp \
	FEManagerT_bridging_1.cpp \
	FEManagerT_bridging_2.cpp \
	FEManagerT_bridging_22.cpp \
	FEManagerT_bridging_3.cpp \
	FEManagerT_bridging_4.cpp \
	GlobalT.cpp \
	IOManager.cpp \
	IOManager_mpi.cpp \
	main.cpp \
	TimeManagerT.cpp \
	CommManagerT.cpp \
	MultiManagerT.cpp

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
	ExecutionManagerT.h \
	FEExecutionManagerT.h \
	FEManagerT.h \
	FEManagerT_mpi.h \
	FEManagerT_bridging.h \
	GlobalT.h \
	IOManager.h \
	IOManager_mpi.h \
	TimeManagerT.h \
	CommManagerT.h \
	MultiManagerT.h

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
