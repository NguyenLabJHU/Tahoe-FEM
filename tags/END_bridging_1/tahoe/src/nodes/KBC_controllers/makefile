# $Id: makefile,v 1.7.6.1 2003-02-10 02:12:23 paklein Exp $
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
	PeriodicNodesT.o \
	KBC_PrescribedT.o

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
	PeriodicNodesT.cpp \
	KBC_PrescribedT.cpp

# headers
HDR = \
	BimaterialK_FieldT.h \
	KBC_ControllerT.h \
	K_FieldT.h \
	MappedPeriodicT.h \
	SymmetricNodesT.h \
	TiedNodesT.h \
	PeriodicNodesT.h \
	KBC_PrescribedT.h

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
