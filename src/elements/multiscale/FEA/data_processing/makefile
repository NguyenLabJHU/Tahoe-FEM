#/* $Id: makefile,v 1.3 2005-08-12 05:38:39 paklein Exp $ */
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	FEA_Data_ProcessorT.o 	\
	FEA_Data_Pro_Diagnostics.o 	

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)

# sources
SRC = \
	FEA_Data_ProcessorT.cpp 	\
	FEA_Data_Pro_Diagnostics.cpp 	

# headers
HDR =	FEA_Data_ProcessorT.h

HDR_LINK = $(HDR:.h=.h_link)

# subdirectories
SUB_DIR = 

subdir_driver: 

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = $(MACRO_DIR)/empty
include $(DEPEND)
