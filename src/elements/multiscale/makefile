# $Id: makefile,v 1.6 2005-08-12 05:38:31 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	StaggeredMultiScaleT.o \
	Traction_and_Body_Force.o \
	FEA_FormatT.o 

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = $(MACRO_DIR)/empty

# sources
SRC = \
	StaggeredMultiScaleT.cpp \
	Traction_and_Body_Force.cpp \
	FEA_FormatT.cpp 

# headers
HDR = \
	StaggeredMultiScaleT.h \
	FEA_FormatT.h 

HDR_LINK = $(HDR:.h=.h_link)

# subdirectories
SUB_DIR = \
	FEA \
	materials \
	continuum

subdir_driver: subdir_loop

# skip missing directory
#SUB_DIR =
#subdir_driver:

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = $(MACRO_DIR)/empty
include $(DEPEND)

