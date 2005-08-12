# $Id: makefile,v 1.4 2005-08-12 05:27:13 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	SubMtx.basics.o \
	SubMtx.init.o \
	SubMtx.initRandom.o \
	SubMtx.instance.o \
	SubMtx.IO.o \
	SubMtx.scalevec.o \
	SubMtx.solve.o \
	SubMtx.solveH.o \
	SubMtx.solveT.o \
	SubMtx.solveupd.o \
	SubMtx.solveupdH.o \
	SubMtx.solveupdT.o \
	SubMtx.sort.o \
	SubMtx.util.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = $(MACRO_DIR)/empty

# sources
SRC = \
	SubMtx.basics.c \
	SubMtx.init.c \
	SubMtx.initRandom.c \
	SubMtx.instance.c \
	SubMtx.IO.c \
	SubMtx.scalevec.c \
	SubMtx.solve.c \
	SubMtx.solveH.c \
	SubMtx.solveT.c \
	SubMtx.solveupd.c \
	SubMtx.solveupdH.c \
	SubMtx.solveupdT.c \
	SubMtx.sort.c \
	SubMtx.util.c

# headers
HDR = \
	SubMtx.h
HDR_LINK = $(HDR:.h=.h_link)

# subdirectories
SUB_DIR = 
#subdir_driver: subdir_loop
# SUB_DIR is empty
subdir_driver: 

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = $(MACRO_DIR)/empty
include $(DEPEND)
