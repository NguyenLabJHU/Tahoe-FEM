# $Id: makefile,v 1.5 2002-04-10 16:48:24 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	BimaterialK_FieldT.o \
	KBC_ControllerT.o \
	K_FieldT.o \
	MappedPeriodicT.o \
	TiedNodesT.o

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	BimaterialK_FieldT.cpp \
	KBC_ControllerT.cpp \
	K_FieldT.cpp \
	MappedPeriodicT.cpp \
	TiedNodesT.cpp

# headers
HDR = \
	BimaterialK_FieldT.h \
	KBC_ControllerT.h \
	K_FieldT.h \
	MappedPeriodicT.h \
	TiedNodesT.h

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
