# $Id: makefile,v 1.3 2001-02-22 23:42:37 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	ZV.basics.o \
	ZV.init.o \
	ZV.instance.o \
	ZV.IO.o \
	ZV.profile.o \
	ZV.util.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	ZV.basics.c \
	ZV.init.c \
	ZV.instance.c \
	ZV.IO.c \
	ZV.profile.c \
	ZV.util.c

# headers
HDR = \
	ZV.h
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
