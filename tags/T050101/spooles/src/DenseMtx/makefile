# $Id: makefile,v 1.3 2001-02-22 23:42:18 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	DenseMtx.basics.o \
	DenseMtx.initAsSubmtx.o \
	DenseMtx.init.o \
	DenseMtx.instance.o \
	DenseMtx.IO.o \
	DenseMtx.permute.o \
	DenseMtx.scale.o \
	DenseMtx.util.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	DenseMtx.basics.c \
	DenseMtx.initAsSubmtx.c \
	DenseMtx.init.c \
	DenseMtx.instance.c \
	DenseMtx.IO.c \
	DenseMtx.permute.c \
	DenseMtx.scale.c \
	DenseMtx.util.c

# headers
HDR = \
	DenseMtx.h
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
