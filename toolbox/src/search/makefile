# $Id: makefile,v 1.5 2001-02-22 23:45:27 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	CartesianGridT.o \
	iGridManagerT.o \
	iGridManager2DT.o \
	iGridManager3DT.o \
	iNodeT.o \
	iPeriodicGrid3DT.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	CartesianGridT.cpp \
	iGridManagerT.cpp \
	iGridManager2DT.cpp \
	iGridManager3DT.cpp \
	iNodeT.cpp \
	iPeriodicGrid3DT.cpp

# headers
HDR = \
	CartesianGridT.h \
	GridManager2DT.h \
	GridManager3DT.h \
	iGridManagerT.h \
	iGridManager2DT.h \
	iGridManager3DT.h \
	iNodeT.h \
	iPeriodicGrid3DT.h
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
