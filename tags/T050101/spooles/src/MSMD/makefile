# $Id: makefile,v 1.3 2001-02-22 23:42:29 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	MSMD.basics.o \
	MSMD.cleanReachSet.o \
	MSMD.eliminate.o \
	MSMD.fillPerms.o \
	MSMD.findInodes.o \
	MSMD.frontETree.o \
	MSMD.init.o \
	MSMD.makeSchurComplement.o \
	MSMD.MSMDinfo.o \
	MSMD.MSMDvtx.o \
	MSMD.order.o \
	MSMD.update.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	MSMD.basics.c \
	MSMD.cleanReachSet.c \
	MSMD.eliminate.c \
	MSMD.fillPerms.c \
	MSMD.findInodes.c \
	MSMD.frontETree.c \
	MSMD.init.c \
	MSMD.makeSchurComplement.c \
	MSMD.MSMDinfo.c \
	MSMD.MSMDvtx.c \
	MSMD.order.c \
	MSMD.update.c

# headers
HDR = \
	MSMD.h
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
