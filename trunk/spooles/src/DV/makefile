# $Id: makefile,v 1.3 2001-02-22 23:42:17 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	DV.basics.o \
	DV.init.o \
	DV.instance.o \
	DV.IO.o \
	DV.profile.o \
	DV.util.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	DV.basics.c \
	DV.init.c \
	DV.instance.c \
	DV.IO.c \
	DV.profile.c \
	DV.util.c

# headers
HDR = \
	DV.h
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
