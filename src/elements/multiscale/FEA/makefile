#/* $Id: makefile,v 1.1 2002-11-22 21:35:39 creigh Exp $ */
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = all_depend.o 

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = all_depend.cpp 

# headers
HDR = FEA.h 
HDR_LINK = $(HDR:.h=.h_link)

# subdirectories
SUB_DIR = \
	matrix \
	integration \
	materials \
	data_processing \
	shape_functions \
	multiscale_elements 

subdir_driver: subdir_loop

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies (none this directory)
DEPEND = /dev/null
include $(DEPEND)
