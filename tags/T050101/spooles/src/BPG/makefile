# $Id: makefile,v 1.3 2001-02-22 23:42:12 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	BPG.basics.o \
	BPG.DM.o \
	BPG.init.o \
	BPG.IO.o \
	BPG.makeGraphs.o \
	BPG.maxFlow.o \
	BPG.pseudo.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	BPG.basics.c \
	BPG.DM.c \
	BPG.init.c \
	BPG.IO.c \
	BPG.makeGraphs.c \
	BPG.maxFlow.c \
	BPG.pseudo.c

# headers
HDR = \
	BPG.h
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
