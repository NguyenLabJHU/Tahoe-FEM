#/* $Id: makefile,v 1.1 2002-11-22 21:38:47 creigh Exp $ */
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = VMS_BCJT.o 		

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)

# sources
SRC = VMS_BCJT.cpp 		

# headers
HDR = \
	FineScaleT.h 	\
	VMS_BCJT.h 		
HDR_LINK = $(HDR:.h=.h_link)

# subdirectories
SUB_DIR = 

subdir_driver: 

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = /dev/null
include $(DEPEND)
