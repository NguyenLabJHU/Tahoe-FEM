# $Id: makefile,v 1.4 2005-08-12 05:27:09 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	A2.basics.o \
	A2.copyEntriesToVector.o \
	A2.init.o \
	A2.instance.o \
	A2.IO.o \
	A2.makeStaircase.o \
	A2.norms.o \
	A2.QRreduce.o \
	A2.sort.o \
	A2.util.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = $(MACRO_DIR)/empty

# sources
SRC = \
	A2.basics.c \
	A2.copyEntriesToVector.c \
	A2.init.c \
	A2.instance.c \
	A2.IO.c \
	A2.makeStaircase.c \
	A2.norms.c \
	A2.QRreduce.c \
	A2.sort.c \
	A2.util.c

# headers
HDR = \
	A2.h
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
