# $Id: makefile,v 1.3 2001-02-22 23:42:24 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	ILUMtx.basics.o \
	ILUMtx.factor.o \
	ILUMtx.init.o \
	ILUMtx.IO.o \
	ILUMtx.misc.o \
	ILUMtx.solve.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	ILUMtx.basics.c \
	ILUMtx.factor.c \
	ILUMtx.init.c \
	ILUMtx.IO.c \
	ILUMtx.misc.c \
	ILUMtx.solve.c

# headers
HDR = \
	ILUMtx.h
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
