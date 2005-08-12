# $Id: makefile,v 1.4 2005-08-12 05:27:12 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	Network.addArc.o \
	Network.augmentPath.o \
	Network.basics.o \
	Network.findAugmentingPath.o \
	Network.findMaxFlow.o \
	Network.findMincut.o \
	Network.init.o \
	Network.IO.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = $(MACRO_DIR)/empty

# sources
SRC = \
	Network.addArc.c \
	Network.augmentPath.c \
	Network.basics.c \
	Network.findAugmentingPath.c \
	Network.findMaxFlow.c \
	Network.findMincut.c \
	Network.init.c \
	Network.IO.c

# headers
HDR = \
	Network.h
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
