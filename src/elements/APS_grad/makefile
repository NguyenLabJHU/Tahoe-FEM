# $Id: makefile,v 1.4 2003-09-03 00:47:26 raregue Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = APS_AssemblyT.o FEA_FormatT.o \
	APS_Traction_and_Body_Force.o

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = APS_AssemblyT.cpp FEA_FormatT.cpp \
	APS_Traction_and_Body_Force.cpp

# headers
HDR = APS_AssemblyT.h  FEA_FormatT.h
HDR_LINK = $(HDR:.h=.h_link)

# subdirectories
SUB_DIR = FEA \
	materials 

subdir_driver: subdir_loop

# skip missing directory
#SUB_DIR =
#subdir_driver:

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = /dev/null
include $(DEPEND)

