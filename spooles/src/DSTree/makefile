# $Id: makefile,v 1.4 2005-08-12 05:27:10 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	DSTree.basics.o \
	DSTree.init.o \
	DSTree.instance.o \
	DSTree.IO.o \
	DSTree.stages.o \
	DSTree.util.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = $(MACRO_DIR)/empty

# sources
SRC = \
	DSTree.basics.c \
	DSTree.init.c \
	DSTree.instance.c \
	DSTree.IO.c \
	DSTree.stages.c \
	DSTree.util.c

# headers
HDR = \
	DSTree.h
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
