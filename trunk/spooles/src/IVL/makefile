# $Id: makefile,v 1.4 2005-08-12 05:27:11 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	IVL.basics.o \
	IVL.init.o \
	IVL.instance.o \
	IVL.IO.o \
	IVL.listmanip.o \
	IVL.misc.o \
	IVL.subIVL.o \
	IVL.util.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = $(MACRO_DIR)/empty

# sources
SRC = \
	IVL.basics.c \
	IVL.init.c \
	IVL.instance.c \
	IVL.IO.c \
	IVL.listmanip.c \
	IVL.misc.c \
	IVL.subIVL.c \
	IVL.util.c

# headers
HDR = \
	IVL.h
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
