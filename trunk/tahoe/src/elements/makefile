# $Id: makefile,v 1.5 2001-08-20 06:45:05 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	ElementT.o \
	ElementBaseT.o \
	ElementListT.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	ElementT.cpp \
	ElementBaseT.cpp \
	ElementListT.cpp

# headers
HDR = \
	ElementT.h \
	ElementBaseT.h \
	ElementListT.h
HDR_LINK = $(HDR:.h=.h_link)

# subdirectories
SUB_DIR = \
	BEM_element \
	cohesive_surface \
	contact \
	continuum \
	lattice \
	rod \
	variable_elements
subdir_driver: subdir_loop
# SUB_DIR is empty
#subdir_driver: 

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = /dev/null
include $(DEPEND)
