# $Id: makefile,v 1.5 2003-01-29 09:10:56 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	MeshFreeElementSupportT.o \
	MeshFreeFractureSupportT.o \
	MeshFreeSSSolidT.o \
	MeshFreeFSSolidT.o \
	D2MeshFreeFSSolidT.o

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	MeshFreeElementSupportT.cpp \
	MeshFreeFractureSupportT.cpp \
	MeshFreeSSSolidT.cpp \
	MeshFreeFSSolidT.cpp \
	D2MeshFreeFSSolidT.cpp

# headers
HDR = \
	MeshFreeElementSupportT.h \
	MeshFreeFractureSupportT.h \
	MeshFreeSSSolidT.h \
	MeshFreeFSSolidT.h \
	D2MeshFreeFSSolidT.h

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
