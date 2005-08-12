#/* $Id: makefile,v 1.3 2005-08-12 05:38:35 paklein Exp $ */
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = all_depend.o 

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = $(MACRO_DIR)/empty

# sources
SRC = all_depend.cpp 

# headers
HDR = FEA.h 
HDR_LINK = $(HDR:.h=.h_link)

# subdirectories
SUB_DIR = \
	matrix \
	integration \
	data_processing \
	shape_functions \
	multiscale_elements 

subdir_driver: subdir_loop

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies (none this directory)
DEPEND = $(MACRO_DIR)/empty
include $(DEPEND)
