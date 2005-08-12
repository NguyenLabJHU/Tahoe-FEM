# $Id: makefile,v 1.4 2005-08-12 05:27:10 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	Chv.assemble.o \
	Chv.basics.o \
	Chv.copy.o \
	Chv.factor.o \
	Chv.findPivot.o \
	Chv.init.o \
	Chv.instance.o \
	Chv.IO.o \
	Chv.search.o \
	Chv.swap.o \
	Chv.update.o \
	Chv.util.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = $(MACRO_DIR)/empty

# sources
SRC = \
	Chv.assemble.c \
	Chv.basics.c \
	Chv.copy.c \
	Chv.factor.c \
	Chv.findPivot.c \
	Chv.init.c \
	Chv.instance.c \
	Chv.IO.c \
	Chv.search.c \
	Chv.swap.c \
	Chv.update.c \
	Chv.util.c

# headers
HDR = \
	Chv.h
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
