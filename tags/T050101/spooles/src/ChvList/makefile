# $Id: makefile,v 1.3 2001-02-22 23:42:14 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	ChvList.basics.o \
	ChvList.init.o \
	ChvList.IO.o \
	ChvList.util.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	ChvList.basics.c \
	ChvList.init.c \
	ChvList.IO.c \
	ChvList.util.c

# headers
HDR = \
	ChvList.h
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
