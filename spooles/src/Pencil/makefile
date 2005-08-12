# $Id: makefile,v 1.4 2005-08-12 05:27:12 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	Pencil.basics.o \
	Pencil.init.o \
	Pencil.IO.o \
	Pencil.mmm.o \
	Pencil.permute.o \
	Pencil.setup.o \
	Pencil.util.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = $(MACRO_DIR)/empty

# sources
SRC = \
	Pencil.basics.c \
	Pencil.init.c \
	Pencil.IO.c \
	Pencil.mmm.c \
	Pencil.permute.c \
	Pencil.setup.c \
	Pencil.util.c

# headers
HDR = \
	Pencil.h
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
