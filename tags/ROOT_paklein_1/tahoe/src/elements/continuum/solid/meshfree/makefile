# $Id: makefile,v 1.4 2001-02-22 23:46:19 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	D2MeshFreeFDElasticT.o \
	MeshFreeElasticT.o \
	MeshFreeElementSupportT.o \
	MeshFreeFDElasticT.o \
	MeshFreeFractureSupportT.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	D2MeshFreeFDElasticT.cpp \
	MeshFreeElasticT.cpp \
	MeshFreeElementSupportT.cpp \
	MeshFreeFDElasticT.cpp \
	MeshFreeFractureSupportT.cpp

# headers
HDR = \
	D2MeshFreeFDElasticT.h \
	MeshFreeElasticT.h \
	MeshFreeElementSupportT.h \
	MeshFreeFDElasticT.h \
	MeshFreeFractureSupportT.h
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
