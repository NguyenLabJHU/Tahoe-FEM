# $Id: makefile,v 1.7 2002-06-27 01:21:05 hspark Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	DiffusionMatListT.o \
	MaterialListT.o \
	MaterialList1DT.o \
	MaterialList2DT.o \
	MaterialList3DT.o \
	StructuralMatListT.o \
	MaterialT.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	DiffusionMatListT.cpp \
	MaterialListT.cpp \
	MaterialList1DT.cpp \
	MaterialList2DT.cpp \
	MaterialList3DT.cpp \
	StructuralMatListT.cpp \
	MaterialT.cpp

# headers
HDR = \
	DiffusionMatListT.h \
	MaterialListT.h \
	MaterialList1DT.h \
	MaterialList2DT.h \
	MaterialList3DT.h \
	StructuralMatListT.h \
	MaterialT.h
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
