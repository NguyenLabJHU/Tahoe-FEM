#/* $Id: makefile,v 1.2 2003-09-22 20:53:14 raregue Exp $ */
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	APS_VariableT.o

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)

# sources
SRC = \
	APS_VariableT.cpp

# headers
HDR =	\
	APS_EnumT.h	\
	APS_VariableT.h	

HDR_LINK = $(HDR:.h=.h_link)

# subdirectories
SUB_DIR = \
	bal_lin \
	plast

subdir_driver: subdir_loop

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = /dev/null
include $(DEPEND)
