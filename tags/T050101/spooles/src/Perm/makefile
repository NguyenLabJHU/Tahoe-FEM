# $Id: makefile,v 1.3 2001-02-22 23:42:31 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	Perm.basics.o \
	Perm.compress.o \
	Perm.init.o \
	Perm.IO.o \
	Perm.util.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	Perm.basics.c \
	Perm.compress.c \
	Perm.init.c \
	Perm.IO.c \
	Perm.util.c

# headers
HDR = \
	Perm.h
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
