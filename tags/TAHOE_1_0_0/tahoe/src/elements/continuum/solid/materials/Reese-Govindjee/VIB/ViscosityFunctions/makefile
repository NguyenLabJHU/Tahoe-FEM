# $Id: makefile,v 1.1 2002-10-04 23:59:03 thao Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	ViscFuncT.o \
	Linear.o \
	QuadraticT.o \
	ConstantT.o \
	BiQuadraticT.o \
	VariViscT.o 

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	ViscFuncT.cpp \
	Linear.cpp \
	QuadraticT.cpp \
	ConstantT.cpp \
	BiQuadraticT.cpp \
	VariViscT.cpp  

# headers
HDR = \
	ViscFuncT.h \
	Linear.h \
	QuadraticT.h \
	ConstantT.h \
	BiQuadraticT.h \
	VariViscT.h 

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
