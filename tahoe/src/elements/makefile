# $Id: makefile,v 1.11 2002-11-26 01:56:37 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	ElementT.o \
	ElementSupportT.o \
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
	ElementSupportT.cpp \
	ElementBaseT.cpp \
	ElementListT.cpp

# headers
HDR = \
	ElementT.h \
	ElementSupportT.h \
	ElementBaseT.h \
	ElementListT.h
HDR_LINK = $(HDR:.h=.h_link)

# subdirectories
SUB_DIR = \
	BEM_element \
	bridging_scale \
	cohesive_surface \
	contact \
	continuum \
	lattice \
	rod \
	adhesion \
	multiscale \
	particle

subdir_driver: subdir_loop
# SUB_DIR is empty
#subdir_driver: 

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = /dev/null
include $(DEPEND)
