# $Id: makefile,v 1.3 2001-02-22 23:42:11 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	BKL.basics.o \
	BKL.evalfcn.o \
	BKL.exhSearch.o \
	BKL.fidmat.o \
	BKL.init.o \
	BKL.util.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	BKL.basics.c \
	BKL.evalfcn.c \
	BKL.exhSearch.c \
	BKL.fidmat.c \
	BKL.init.c \
	BKL.util.c

# headers
HDR = \
	BKL.h
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
