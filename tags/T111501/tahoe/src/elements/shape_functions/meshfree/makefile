# $Id: makefile,v 1.5 2001-06-13 20:53:55 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	MeshFreeShapeFunctionT.o \
	MeshFreeSupportT.o \
	MeshFreeSupport2DT.o \
	MeshFreeSupport3DT.o \
	MeshFreeSurfaceShapeT.o \
	MeshFreeSurfaceSupportT.o \
	MeshFreeT.o \
	MLSSolverT.o \
	OrthoMLSSolverT.o \
	OrthoMLS2DT.o \
	OrthoMLS3DT.o \
	SamplingSurfaceT.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	MeshFreeShapeFunctionT.cpp \
	MeshFreeSupportT.cpp \
	MeshFreeSupport2DT.cpp \
	MeshFreeSupport3DT.cpp \
	MeshFreeSurfaceShapeT.cpp \
	MeshFreeSurfaceSupportT.cpp \
	MeshFreeT.cpp \
	MLSSolverT.cpp \
	OrthoMLSSolverT.cpp \
	OrthoMLS2DT.cpp \
	OrthoMLS3DT.cpp \
	SamplingSurfaceT.cpp

# headers
HDR = \
	MeshFreeShapeFunctionT.h \
	MeshFreeSupportT.h \
	MeshFreeSupport2DT.h \
	MeshFreeSupport3DT.h \
	MeshFreeSurfaceShapeT.h \
	MeshFreeSurfaceSupportT.h \
	MeshFreeT.h \
	MLSSolverT.h \
	OrthoMLSSolverT.h \
	OrthoMLS2DT.h \
	OrthoMLS3DT.h \
	SamplingSurfaceT.h
HDR_LINK = $(HDR:.h=.h_link)

# subdirectories
SUB_DIR = \
	basis \
	D2 \
	windows
subdir_driver: subdir_loop
# SUB_DIR is empty
#subdir_driver: 

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = /dev/null
include $(DEPEND)
