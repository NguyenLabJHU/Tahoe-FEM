# $Id: makefile,v 1.4 2001-02-22 23:46:20 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	RodT.o \
	UnConnectedRodT.o \
	VirtualRodT.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	RodT.cpp \
	UnConnectedRodT.cpp \
	VirtualRodT.cpp

# headers
HDR = \
	RodT.h \
	UnConnectedRodT.h \
	VirtualRodT.h
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
