# $Id: makefile,v 1.10 2002-09-22 19:54:59 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	VTKUGridT.o \
	VTKBodyT.o \
	VTKBodyDataT.o \
	VTKConsoleT.o \
	VTKFrameT.o \
	VTKMappedIdFilterT.o \
	vtkTahoeGlyph3D.o \
	main.o \
	gl2ps.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	VTKBodyT.cpp \
	VTKUGridT.cpp \
	VTKBodyDataT.cpp \
	VTKConsoleT.cpp \
	VTKFrameT.cpp \
	VTKMappedIdFilterT.cpp \
	vtkTahoeGlyph3D.cpp \
	main.cpp \
	gl2ps.c

# headers
HDR = \
	VTKBodyT.h \
	VTKUGridT.h \
	VTKBodyDataT.h \
	VTKConsoleT.h \
	VTKFrameT.h \
	VTKMappedIdFilterT.h \
	vtkTahoeGlyph3D.h \
	gl2ps.h

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
