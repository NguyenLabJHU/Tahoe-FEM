# $Id: makefile,v 1.3 2001-02-22 23:42:23 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	I2Ohash.basics.o \
	I2Ohash.init.o \
	I2Ohash.IO.o \
	I2Ohash.util.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	I2Ohash.basics.c \
	I2Ohash.init.c \
	I2Ohash.IO.c \
	I2Ohash.util.c

# headers
HDR = \
	I2Ohash.h
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
