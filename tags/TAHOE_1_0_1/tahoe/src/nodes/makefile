# $Id: makefile,v 1.8 2002-06-08 20:20:48 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	BasicFieldT.o \
	FieldT.o \
	FieldSupportT.o \
	NodeManagerT.o

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	BasicFieldT.cpp \
	FieldT.cpp \
	FieldSupportT.cpp \
	NodeManagerT.cpp

# headers
HDR = \
	BasicFieldT.h \
	FieldT.h \
	FieldSupportT.h \
	NodeManagerT.h

HDR_LINK = $(HDR:.h=.h_link)

# subdirectories
SUB_DIR = \
	FBC_controllers \
	KBC_controllers \
	x_dof
subdir_driver: subdir_loop
# SUB_DIR is empty
#subdir_driver: 

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = /dev/null
include $(DEPEND)
