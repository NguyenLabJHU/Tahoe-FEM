# $Id: makefile,v 1.5 2001-02-22 23:45:16 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	BandedLAdMatrixT.o \
	dMatrixEXT.o \
	dMatrixT.o \
	dSymMatrixT.o \
	LAdMatrixT.o \
	TriDiagdMatrixT.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	BandedLAdMatrixT.cpp \
	dMatrixEXT.cpp \
	dMatrixT.cpp \
	dSymMatrixT.cpp \
	LAdMatrixT.cpp \
	TriDiagdMatrixT.cpp

# headers
HDR = \
	BandedLAdMatrixT.h \
	dMatrixEXT.h \
	dMatrixT.h \
	dSymMatrixT.h \
	LAdMatrixT.h \
	nMatrixT.h \
	TriDiagdMatrixT.h
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
