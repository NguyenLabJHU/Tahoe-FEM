# $Id: makefile,v 1.3 2001-02-22 23:42:19 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	EGraph.basics.o \
	EGraph.init.o \
	EGraph.IO.o \
	EGraph.misc.o \
	EGraph.mkAdjGraph.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	EGraph.basics.c \
	EGraph.init.c \
	EGraph.IO.c \
	EGraph.misc.c \
	EGraph.mkAdjGraph.c

# headers
HDR = \
	EGraph.h
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
