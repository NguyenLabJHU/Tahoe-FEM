#/* $Id: makefile,v 1.2 2002-11-23 16:32:41 paklein Exp $ */
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	VMF_MaterialT.o 	

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)

# sources
SRC = \
	VMF_MaterialT.cpp 	

# headers
HDR =	VMF_MaterialT.h	\
			E_Pr_MatlT.h \
			Iso_MatlT.h \
			BCJ_MatlT.h

HDR_LINK = $(HDR:.h=.h_link)

# subdirectories
SUB_DIR = 

subdir_driver: 

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = /dev/null
include $(DEPEND)
