# $Id: makefile,v 1.3 2001-02-22 23:42:32 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	SemiImplMtx.basics.o \
	SemiImplMtx.init.o \
	SemiImplMtx.IO.o \
	SemiImplMtx.solve.o \
	SemiImplMtx.util.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	SemiImplMtx.basics.c \
	SemiImplMtx.init.c \
	SemiImplMtx.IO.c \
	SemiImplMtx.solve.c \
	SemiImplMtx.util.c

# headers
HDR = \
	SemiImplMtx.h
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
