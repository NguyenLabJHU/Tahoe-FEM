# $Id: makefile,v 1.3 2001-02-22 23:42:25 paklein Exp $
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
#DEP = /dev/null

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
DEPEND = /dev/null
include $(DEPEND)
