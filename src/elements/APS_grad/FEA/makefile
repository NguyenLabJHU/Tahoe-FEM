#/* $Id: makefile,v 1.1 2003-07-10 17:18:01 raregue Exp $ */
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = 

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = 

# headers
HDR = FEA.h 
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
