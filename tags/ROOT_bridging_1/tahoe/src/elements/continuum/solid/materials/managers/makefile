# $Id: makefile,v 1.9 2003-01-30 00:43:39 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	SolidMatListT.o \
	SolidMatList1DT.o \
	SolidMatList2DT.o \
	SolidMatList3DT.o \
	SolidT.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	SolidMatListT.cpp \
	SolidMatList1DT.cpp \
	SolidMatList2DT.cpp \
	SolidMatList3DT.cpp \
	SolidT.cpp

# headers
HDR = \
	SolidMatListT.h \
	SolidMatList1DT.h \
	SolidMatList2DT.h \
	SolidMatList3DT.h \
	SolidT.h

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
