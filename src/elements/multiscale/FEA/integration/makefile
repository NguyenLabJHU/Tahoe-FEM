#/* $Id: makefile,v 1.5 2005-08-12 05:38:44 paklein Exp $ */
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	FEA_IntegrationT.o \
	FEA_SurfIntegrationT.o

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)

# sources
SRC = \
	FEA_IntegrationT.cpp \
	FEA_SurfIntegrationT.cpp

# headers
HDR = \
	FEA_IntegrationT.h \
	FEA_SurfIntegrationT.h
		 
HDR_LINK = $(HDR:.h=.h_link)

# subdirectories
SUB_DIR = 

subdir_driver: 

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = $(MACRO_DIR)/empty
include $(DEPEND)
