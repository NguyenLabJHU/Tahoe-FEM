# $Id: makefile,v 1.6 2001-06-23 06:25:44 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	GaussianWindowT.o \
	RectCubicSplineWindowT.o \
	CubicSplineWindowT.o \
	RectGaussianWindowT.o

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	GaussianWindowT.cpp \
	RectCubicSplineWindowT.cpp \
	CubicSplineWindowT.cpp \
	RectGaussianWindowT.cpp

# headers
HDR = \
	GaussianWindowT.h \
	RectCubicSplineWindowT.h \
	RectGaussianWindowT.h \
	CubicSplineWindowT.h \
	WindowT.h

HDR_LINK = $(HDR:.h=.h_link)

# subdirectories
#SUB_DIR =
#subdir_driver: subdir_loop

# SUB_DIR is empty
subdir_driver: 

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = /dev/null
include $(DEPEND)
