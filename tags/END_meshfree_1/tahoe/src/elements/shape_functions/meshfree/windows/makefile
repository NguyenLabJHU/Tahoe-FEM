# $Id: makefile,v 1.3 2001-06-19 23:25:05 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	GaussianWindowT.o \
	RectangularGaussianWindowT.o

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	GaussianWindowT.cpp \
	RectangularGaussianWindowT.cpp

# headers
HDR = \
	GaussianWindowT.h \
	RectangularGaussianWindowT.h \
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
