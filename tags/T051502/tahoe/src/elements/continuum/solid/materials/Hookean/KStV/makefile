# $Id: makefile,v 1.5 2001-07-03 01:35:09 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	FDKStV.o \
	SSKStV.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	FDKStV.cpp \
	SSKStV.cpp

# headers
HDR = \
	FDKStV.h \
	SSKStV.h
HDR_LINK = $(HDR:.h=.h_link)

# subdirectories
SUB_DIR = \
	2D
subdir_driver: subdir_loop
# SUB_DIR is empty
#subdir_driver: 

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = /dev/null
include $(DEPEND)
