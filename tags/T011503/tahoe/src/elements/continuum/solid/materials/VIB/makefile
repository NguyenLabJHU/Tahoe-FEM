# $Id: makefile,v 1.4 2001-02-22 23:46:33 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	VIB.o \
	VIB_E_MatT.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	VIB.cpp \
	VIB_E_MatT.cpp

# headers
HDR = \
	VIB_E_MatT.h \
	VIB.h
HDR_LINK = $(HDR:.h=.h_link)

# subdirectories
SUB_DIR = \
	circlepoints \
	D2VIB \
	sphericalpoints \
	2D \
	3D
subdir_driver: subdir_loop
# SUB_DIR is empty
#subdir_driver: 

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = /dev/null
include $(DEPEND)
