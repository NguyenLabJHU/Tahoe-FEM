# $Id: makefile,v 1.1 2001-07-18 21:29:45 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	ABAQUS_VUMAT_BCJ.o \
	ABAQUS_VUMAT_BaseT.o \
	vumatadiab.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	ABAQUS_VUMAT_BCJ.cpp \
	ABAQUS_VUMAT_BaseT.cpp \
	vumatadiab.c

# headers
HDR = \
	ABAQUS_VUMAT_BCJ.h \
	ABAQUS_VUMAT_BaseT.h
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
