#/* $Id: makefile,v 1.7 2005-08-12 05:39:01 paklein Exp $ */
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	VMS_BCJT.o \
	VMS_BCJ_XT.o \
	VMS_EZ3T.o	\
	VMS_EZT.o

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)

# sources
SRC = \
	VMS_BCJT.cpp	\
	VMS_BCJ_XT.cpp	\
	VMS_EZ3T.cpp	\
	VMS_EZT.cpp

# headers
HDR = \
	FineScaleT.h 	\
	VMS_BCJT.h 	\
	VMS_BCJ_XT.h 	\
	VMS_EZ3T.h	\
	VMS_EZT.h

HDR_LINK = $(HDR:.h=.h_link)

# subdirectories
SUB_DIR = 

subdir_driver: 

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = $(MACRO_DIR)/empty
include $(DEPEND)
