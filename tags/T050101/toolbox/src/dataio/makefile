# $Id: makefile,v 1.5 2001-02-22 23:45:20 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	IOBaseT.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	IOBaseT.cpp

# headers
HDR = \
	IOBaseT.h
HDR_LINK = $(HDR:.h=.h_link)

# subdirectories
SUB_DIR = \
	database \
	input \
	output
subdir_driver: subdir_loop
# SUB_DIR is empty
#subdir_driver: 

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = /dev/null
include $(DEPEND)
