# $Id: makefile,v 1.3 2001-02-22 23:42:22 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	Graph.basics.o \
	Graph.compress.o \
	Graph.equivMap.o \
	Graph.expand.o \
	Graph.fillFromOffsets.o \
	Graph.init.o \
	Graph.IO.o \
	Graph.setListsFromOffsets.o \
	Graph.util.o \
	Graph.wirebasket.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	Graph.basics.c \
	Graph.compress.c \
	Graph.equivMap.c \
	Graph.expand.c \
	Graph.fillFromOffsets.c \
	Graph.init.c \
	Graph.IO.c \
	Graph.setListsFromOffsets.c \
	Graph.util.c \
	Graph.wirebasket.c

# headers
HDR = \
	Graph.h
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
