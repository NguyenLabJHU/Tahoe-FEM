#/* $Id: makefile,v 1.2 2003-11-21 22:54:28 paklein Exp $ */
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	APS_FEA_Data_ProcessorT.o \
	APS_FEA_Data_Processor_SurfT.o

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)

# sources
SRC = \
	APS_FEA_Data_ProcessorT.cpp \
	APS_FEA_Data_Processor_SurfT.cpp	

# headers
HDR = \
	APS_FEA_Data_ProcessorT.h \
	APS_FEA_Data_Processor_SurfT.h

HDR_LINK = $(HDR:.h=.h_link)

# subdirectories
SUB_DIR = 

subdir_driver: 

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = /dev/null
include $(DEPEND)
