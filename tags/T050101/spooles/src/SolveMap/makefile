# $Id: makefile,v 1.3 2001-02-22 23:42:33 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	SolveMap.basics.o \
	SolveMap.init.o \
	SolveMap.instance.o \
	SolveMap.IO.o \
	SolveMap.maps.o \
	SolveMap.setup.o \
	SolveMap.util.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	SolveMap.basics.c \
	SolveMap.init.c \
	SolveMap.instance.c \
	SolveMap.IO.c \
	SolveMap.maps.c \
	SolveMap.setup.c \
	SolveMap.util.c

# headers
HDR = \
	SolveMap.h
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
