# $Id: makefile,v 1.5.2.1 2001-10-09 13:19:57 sawimme Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	DuNodeManager.o \
	DynNodeManager.o \
	FDDynNodeManagerT.o \
	FDNodeManager.o \
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
	NodeManagerPrimitive.cpp \
	NodeManagerT.cpp

# headers
HDR = \
	DuNodeManager.h \
	DynNodeManager.h \
	FDDynNodeManagerT.h \
	FDNodeManager.h \
	NodeManagerPrimitive.h \
	NodeManagerT.h
HDR_LINK = $(HDR:.h=.h_link)

# subdirectories
SUB_DIR = \
	FBC_controllers \
	KBC_controllers \
	x_dof
subdir_driver: subdir_loop
# SUB_DIR is empty
#subdir_driver: 

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = /dev/null
include $(DEPEND)
