# $Id: makefile,v 1.6 2003-02-22 23:53:44 paklein Exp $
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
	TriDiagdMatrixT.o \
	dSPMatrixT.o

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
	TriDiagdMatrixT.cpp \
	dSPMatrixT.cpp

# headers
HDR = \
	BandedLAdMatrixT.h \
	dMatrixEXT.h \
	dMatrixT.h \
	dSymMatrixT.h \
	LAdMatrixT.h \
	nMatrixT.h \
	TriDiagdMatrixT.h \
	dSPMatrixT.h

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
