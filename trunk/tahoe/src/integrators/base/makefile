# $Id: makefile,v 1.4 2001-02-22 23:46:09 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	ControllerT.o \
	eControllerT.o \
	nControllerT.o \
	nDDtControllerT.o \
	nDtControllerT.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	ControllerT.cpp \
	eControllerT.cpp \
	nControllerT.cpp \
	nDDtControllerT.cpp \
	nDtControllerT.cpp

# headers
HDR = \
	ControllerT.h \
	eControllerT.h \
	nControllerT.h \
	nDDtControllerT.h \
	nDtControllerT.h
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
