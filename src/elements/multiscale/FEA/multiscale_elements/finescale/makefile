#/* $Id: makefile,v 1.2 2002-12-12 21:42:48 creigh Exp $ */
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	VMS_BCJT.o \
	VMS_EZT.o

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)

# sources
SRC = \
	VMS_BCJT.cpp	\
	VMS_EZT.cpp		

# headers
HDR = \
	FineScaleT.h 	\
	VMS_BCJT.h 	\
  VMS_EZT.h	
HDR_LINK = $(HDR:.h=.h_link)

# subdirectories
SUB_DIR = 

subdir_driver: 

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = /dev/null
include $(DEPEND)
