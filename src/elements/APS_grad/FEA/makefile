#/* $Id: makefile,v 1.4 2003-09-22 20:53:12 raregue Exp $ */
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ =  all_depend.o

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC =  all_depend.cpp

# headers
HDR = APS_FEA.h
HDR_LINK = $(HDR:.h=.h_link)

# subdirectories
SUB_DIR = \
	data_process \
	coupled_eqs 

subdir_driver: subdir_loop

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies (none this directory)
DEPEND = /dev/null
include $(DEPEND)

