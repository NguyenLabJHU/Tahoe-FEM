# $Id: makefile,v 1.4 2001-02-22 23:46:52 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	CCNSMatrixT.o \
	CCSMatrixT.o \
	DiagonalMatrixT.o \
	FullMatrixT.o \
	GlobalMatrixT.o \
	SLUMatrix.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	CCNSMatrixT.cpp \
	CCSMatrixT.cpp \
	DiagonalMatrixT.cpp \
	FullMatrixT.cpp \
	GlobalMatrixT.cpp \
	SLUMatrix.cpp

# headers
HDR = \
	CCNSMatrixT.h \
	CCSMatrixT.h \
	DiagonalMatrixT.h \
	FullMatrixT.h \
	GlobalMatrixT.h \
	SLUMatrix.h \
	superlu.h
HDR_LINK = $(HDR:.h=.h_link)

# subdirectories
SUB_DIR = \
	aztec \
	SPOOLES
subdir_driver: subdir_loop
# SUB_DIR is empty
#subdir_driver: 

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = /dev/null
include $(DEPEND)
