#/* $Id: makefile,v 1.1 2002-11-22 21:39:32 creigh Exp $ */
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = VMF_Virtual_Work_EqT.o 

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)

# sources
SRC = VMF_Virtual_Work_EqT.cpp 

# headers
HDR = VMF_Virtual_Work_EqT.h	\
			CoarseScaleT.h 
HDR_LINK = $(HDR:.h=.h_link)

# subdirectories
SUB_DIR = 

subdir_driver: 

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = /dev/null
include $(DEPEND)
