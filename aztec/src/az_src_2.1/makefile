# $Id: makefile,v 1.1 2002-03-25 17:46:17 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
include object_list

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
include source_list

# headers
HDR = \
	az_aztec_defs.h \
	az_aztec.h
HDR_LINK = $(HDR:.h=.h_link)

# subdirectories
SUB_DIR = 
#subdir_driver: subdir_loop
# SUB_DIR is empty
subdir_driver: 

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = /dev/null
include $(DEPEND)
