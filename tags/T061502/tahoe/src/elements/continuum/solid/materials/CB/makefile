# $Id: makefile,v 1.4 2001-02-22 23:46:24 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	BondLatticeT.o \
	CBLatticeT.o \
	LJFCC111.o \
	LJTr2D.o \
	SWDiamond100.o \
	SWDiamond110.o \
	SWMaterial2D.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	BondLatticeT.cpp \
	CBLatticeT.cpp \
	LJFCC111.cpp \
	LJTr2D.cpp \
	SWDiamond100.cpp \
	SWDiamond110.cpp \
	SWMaterial2D.cpp

# headers
HDR = \
	BondLatticeT.h \
	CBLatticeT.h \
	LJFCC111.h \
	LJTr2D.h \
	SWDiamond100.h \
	SWDiamond110.h \
	SWMaterial2D.h
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
