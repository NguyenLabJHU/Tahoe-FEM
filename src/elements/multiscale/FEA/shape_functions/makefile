#/* $Id: makefile,v 1.1 2002-11-22 21:45:32 creigh Exp $ */
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = all_depend_2.o

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)

# sources
SRC = all_depend_2.cpp

# headers
HDR =	FEA_ShapeFunctionT.h 	
	HDR_LINK = $(HDR:.h=.h_link)

# subdirectories
SUB_DIR = 

subdir_driver: 

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies - none this directory (no object code)
DEPEND = /dev/null
include $(DEPEND)
