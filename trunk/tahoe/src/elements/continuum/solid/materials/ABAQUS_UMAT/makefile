# $Id: makefile,v 1.4 2001-02-22 23:46:23 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	ABAQUS_BCJ.o \
	ABAQUS_UMAT_BaseT.o \
	cycdmg7.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	ABAQUS_BCJ.cpp \
	ABAQUS_UMAT_BaseT.cpp \
	cycdmg7.c

# headers
HDR = \
	ABAQUS_BCJ.h \
	ABAQUS_UMAT_BaseT.h
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
