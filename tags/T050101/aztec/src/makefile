# $Id: makefile,v 1.4 2001-02-22 23:40:03 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
#OBJ = 
#OBJ_LINK = $(OBJ)
#DEP = $(OBJ:.o=.d)
OBJ = dummy
OBJ_LINK =
DEP = /dev/null

# sources
#SRC = 

# headers
HDR = \
	Environment.h
HDR_LINK = $(HDR:.h=.h_link)

# subdirectories
SUB_DIR = \
	az_src \
	extras \
	LINPACK_f2c \
	Y12M_f2c
subdir_driver: subdir_loop
# SUB_DIR is empty
#subdir_driver: 

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = /dev/null
include $(DEPEND)
