# $Id: makefile,v 1.5 2001-06-22 22:33:51 hspark Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	GaussianWindowT.o \
        RectangularCubicSplineWindowT.o \
        CubicSplineWindowT.o \
	RectangularGaussianWindowT.o

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	GaussianWindowT.cpp \
	RectangularCubicSplineWindowT.cpp \
        CubicSplineWindowT.cpp \
	RectangularGaussianWindowT.cpp

# headers
HDR = \
	GaussianWindowT.h \
	RectangularCubicSplineWindowT.h \
	RectangularGaussianWindowT.h \
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
