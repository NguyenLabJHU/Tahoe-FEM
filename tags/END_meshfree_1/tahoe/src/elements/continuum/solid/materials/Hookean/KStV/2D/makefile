# $Id: makefile,v 1.4 2001-02-22 23:46:28 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	FDKStV2D.o \
	KStV2D.o \
	SSKStV2D.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	FDKStV2D.cpp \
	KStV2D.cpp \
	SSKStV2D.cpp

# headers
HDR = \
	FDKStV2D.h \
	KStV2D.h \
	SSKStV2D.h
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
