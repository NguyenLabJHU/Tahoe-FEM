#/* $Id: makefile,v 1.1 2003-07-10 17:21:56 raregue Exp $ */
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = APS_Bal_EqT.o 

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)

# sources
SRC = APS_Bal_EqT.cpp 

# headers
HDR = APS_Bal_EqT.h	\
	BalLinMomT.h 

HDR_LINK = $(HDR:.h=.h_link)

# subdirectories
SUB_DIR = 

subdir_driver: 

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = /dev/null
include $(DEPEND)
