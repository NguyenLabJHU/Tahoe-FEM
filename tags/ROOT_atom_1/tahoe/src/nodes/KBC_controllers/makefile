# $Id: makefile,v 1.7 2002-10-20 22:49:29 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	BimaterialK_FieldT.o \
	KBC_ControllerT.o \
	K_FieldT.o \
	MappedPeriodicT.o \
	SymmetricNodesT.o \
	TiedNodesT.o \
	PeriodicNodesT.o

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
	SymmetricNodesT.cpp \
	TiedNodesT.cpp \
	PeriodicNodesT.cpp

# headers
HDR = \
	BimaterialK_FieldT.h \
	KBC_ControllerT.h \
	K_FieldT.h \
	MappedPeriodicT.h \
	SymmetricNodesT.h \
	TiedNodesT.h \
	PeriodicNodesT.h

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
