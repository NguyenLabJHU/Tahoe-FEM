# $Id: makefile,v 1.3 2001-02-22 23:42:29 paklein Exp $
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
#DEP = /dev/null

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
DEPEND = /dev/null
include $(DEPEND)
