# $Id: makefile,v 1.5 2001-04-27 16:47:58 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	DiffusionMatListT.o \
	MaterialListT.o \
	MaterialList2DT.o \
	MaterialList3DT.o \
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
	MaterialList2DT.cpp \
	MaterialList3DT.cpp \
	MaterialT.cpp

# headers
HDR = \
	DiffusionMatListT.h \
	MaterialListT.h \
	MaterialList2DT.h \
	MaterialList3DT.h \
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
