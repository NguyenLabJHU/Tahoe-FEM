# $Id: makefile,v 1.4 2001-02-22 23:46:45 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	BondLengthsT.o \
	LengthsAndAnglesT.o \
	ModCBSolverT.o \
	ModCB2DT.o \
	ModCB3DT.o \
	PTHT2BodyT.o \
	PTHT3BodyT.o \
	SW2BodyT.o \
	SW3BodyT.o \
	ThreeBodyT.o \
	TwoBodyT.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	BondLengthsT.cpp \
	LengthsAndAnglesT.cpp \
	ModCBSolverT.cpp \
	ModCB2DT.cpp \
	ModCB3DT.cpp \
	PTHT2BodyT.cpp \
	PTHT3BodyT.cpp \
	SW2BodyT.cpp \
	SW3BodyT.cpp \
	ThreeBodyT.cpp \
	TwoBodyT.cpp

# headers
HDR = \
	BondLengthsT.h \
	LengthsAndAnglesT.h \
	ModCBSolverT.h \
	ModCB2DT.h \
	ModCB3DT.h \
	PTHT2BodyT.h \
	PTHT3BodyT.h \
	SW2BodyT.h \
	SW3BodyT.h \
	ThreeBodyT.h \
	TwoBodyT.h
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
