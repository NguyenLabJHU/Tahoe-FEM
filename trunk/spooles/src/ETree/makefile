# $Id: makefile,v 1.4 2005-08-12 05:27:11 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	ETree.basics.o \
	ETree.compress.o \
	ETree.init.o \
	ETree.initFromSubtree.o \
	ETree.instance.o \
	ETree.IO.o \
	ETree.justify.o \
	ETree.maps.o \
	ETree.metrics.o \
	ETree.ms.o \
	ETree.permute.o \
	ETree.semi.o \
	ETree.storage.o \
	ETree.transform.o \
	ETree.util.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = $(MACRO_DIR)/empty

# sources
SRC = \
	ETree.basics.c \
	ETree.compress.c \
	ETree.init.c \
	ETree.initFromSubtree.c \
	ETree.instance.c \
	ETree.IO.c \
	ETree.justify.c \
	ETree.maps.c \
	ETree.metrics.c \
	ETree.ms.c \
	ETree.permute.c \
	ETree.semi.c \
	ETree.storage.c \
	ETree.transform.c \
	ETree.util.c

# headers
HDR = \
	ETree.h
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
