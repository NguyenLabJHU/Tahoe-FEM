# $Id: makefile,v 1.4 2001-02-22 23:46:57 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	D2MeshFreeShapeFunctionT.o \
	D2MeshFreeSupportT.o \
	D2MeshFreeSupport2DT.o \
	D2OrthoMLSSolverT.o \
	D2OrthoMLS1DT.o \
	D2OrthoMLS2DT.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	D2MeshFreeShapeFunctionT.cpp \
	D2MeshFreeSupportT.cpp \
	D2MeshFreeSupport2DT.cpp \
	D2OrthoMLSSolverT.cpp \
	D2OrthoMLS1DT.cpp \
	D2OrthoMLS2DT.cpp

# headers
HDR = \
	D2MeshFreeShapeFunctionT.h \
	D2MeshFreeSupportT.h \
	D2MeshFreeSupport2DT.h \
	D2OrthoMLSSolverT.h \
	D2OrthoMLS1DT.h \
	D2OrthoMLS2DT.h
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
