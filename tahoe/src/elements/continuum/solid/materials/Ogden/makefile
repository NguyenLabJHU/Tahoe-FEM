# $Id: makefile,v 1.5 2001-02-22 23:46:29 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	OgdenIsotropicT.o \
	QuadLogOgden2DT.o \
	QuadLogOgden3DT.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	OgdenIsotropicT.cpp \
	QuadLogOgden2DT.cpp \
	QuadLogOgden3DT.cpp

# headers
HDR = \
	OgdenIsotropicT.h \
	QuadLogOgden2DT.h \
	QuadLogOgden3DT.h
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
