# $Id: makefile,v 1.5 2001-08-20 15:56:58 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	DomainIntegrationT.o \
	ParentDomainT.o \
	ShapeFunctionT.o \
	SimoShapeFunctionT.o \
	SurfaceShapeT.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	DomainIntegrationT.cpp \
	ParentDomainT.cpp \
	ShapeFunctionT.cpp \
	SimoShapeFunctionT.cpp \
	SurfaceShapeT.cpp

# headers
HDR = \
	DomainIntegrationT.h \
	ParentDomainT.h \
	ShapeFunctionT.h \
	SimoShapeFunctionT.h \
	SurfaceShapeT.h
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
