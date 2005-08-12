#/* $Id: makefile,v 1.3 2005-08-12 05:38:52 paklein Exp $ */
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	VMS_VariableT.o

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)

# sources
SRC = \
	VMS_VariableT.cpp

# headers
HDR =	\
	VMS.h	\
	VMS_EnumT.h	\
	VMS_VariableT.h	

HDR_LINK = $(HDR:.h=.h_link)

# subdirectories
SUB_DIR = \
	coarsescale \
	finescale

subdir_driver: subdir_loop

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = $(MACRO_DIR)/empty
include $(DEPEND)
