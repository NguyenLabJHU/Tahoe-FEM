# $Id: makefile,v 1.4 2001-02-22 23:46:30 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	QuadLog2D.o \
	QuadLog3D.o \
	SimoIso2D.o \
	SimoIso3D.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	QuadLog2D.cpp \
	QuadLog3D.cpp \
	SimoIso2D.cpp \
	SimoIso3D.cpp

# headers
HDR = \
	QuadLog2D.h \
	QuadLog3D.h \
	SimoIso2D.h \
	SimoIso3D.h
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
