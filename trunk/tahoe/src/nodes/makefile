# $Id: makefile,v 1.4 2001-02-22 23:46:48 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	DuNodeManager.o \
	DynNodeManager.o \
	FDDynNodeManagerT.o \
	FDNodeManager.o \
	LinStaticNodes.o \
	NodeManagerPrimitive.o \
	NodeManagerT.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	DuNodeManager.cpp \
	DynNodeManager.cpp \
	FDDynNodeManagerT.cpp \
	FDNodeManager.cpp \
	LinStaticNodes.cpp \
	NodeManagerPrimitive.cpp \
	NodeManagerT.cpp

# headers
HDR = \
	DuNodeManager.h \
	DynNodeManager.h \
	FDDynNodeManagerT.h \
	FDNodeManager.h \
	LinStaticNodes.h \
	NodeManagerPrimitive.h \
	NodeManagerT.h
HDR_LINK = $(HDR:.h=.h_link)

# subdirectories
SUB_DIR = \
	FBC_controllers \
	KBC_controllers \
	variable_nodes \
	x_dof
subdir_driver: subdir_loop
# SUB_DIR is empty
#subdir_driver: 

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = /dev/null
include $(DEPEND)
