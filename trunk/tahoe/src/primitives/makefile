# $Id: makefile,v 1.5 2002-01-09 12:00:44 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	ElementCardT.o \
	ElementMatrixT.o \
	FBC_CardT.o \
	IC_CardT.o \
	KBC_CardT.o \
	LoadTime.o \
	SWDataT.o \
	TimeSequence.o \
	Traction_CardT.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	ElementCardT.cpp \
	ElementMatrixT.cpp \
	FBC_CardT.cpp \
	IC_CardT.cpp \
	KBC_CardT.cpp \
	LoadTime.cpp \
	SWDataT.cpp \
	TimeSequence.cpp \
	Traction_CardT.cpp

# headers
HDR = \
	ElementCardT.h \
	ElementMatrixT.h \
	FBC_CardT.h \
	IC_CardT.h \
	KBC_CardT.h \
	LoadTime.h \
	MapSetT.h \
	SWDataT.h \
	TimeSequence.h \
	Traction_CardT.h
HDR_LINK = $(HDR:.h=.h_link)

# subdirectories
SUB_DIR = \
	globalmatrix
subdir_driver: subdir_loop
# SUB_DIR is empty
#subdir_driver: 

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = /dev/null
include $(DEPEND)
