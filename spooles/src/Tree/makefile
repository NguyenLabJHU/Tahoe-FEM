# $Id: makefile,v 1.4 2005-08-12 05:27:13 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	Tree.basics.o \
	Tree.compress.o \
	Tree.draw.o \
	Tree.getCoords.o \
	Tree.init.o \
	Tree.instance.o \
	Tree.IO.o \
	Tree.justify.o \
	Tree.maximizeGain.o \
	Tree.metrics.o \
	Tree.perms.o \
	Tree.permute.o \
	Tree.subtree.o \
	Tree.util.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = $(MACRO_DIR)/empty

# sources
SRC = \
	Tree.basics.c \
	Tree.compress.c \
	Tree.draw.c \
	Tree.getCoords.c \
	Tree.init.c \
	Tree.instance.c \
	Tree.IO.c \
	Tree.justify.c \
	Tree.maximizeGain.c \
	Tree.metrics.c \
	Tree.perms.c \
	Tree.permute.c \
	Tree.subtree.c \
	Tree.util.c

# headers
HDR = \
	Tree.h
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
