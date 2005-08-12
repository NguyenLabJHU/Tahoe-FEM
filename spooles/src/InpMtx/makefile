# $Id: makefile,v 1.4 2005-08-12 05:27:12 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	InpMtx.adjForATA.o \
	InpMtx.basics.o \
	InpMtx.extract.o \
	InpMtx.fullAdj.o \
	InpMtx.gmmm.o \
	InpMtx.gmvm.o \
	InpMtx.init.o \
	InpMtx.input.o \
	InpMtx.instance.o \
	InpMtx.IO.o \
	InpMtx.map.o \
	InpMtx.mvm.o \
	InpMtx.mvmVector.o \
	InpMtx.permute.o \
	InpMtx.profile.o \
	InpMtx.support.o \
	InpMtx.util.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = $(MACRO_DIR)/empty

# sources
SRC = \
	InpMtx.adjForATA.c \
	InpMtx.basics.c \
	InpMtx.extract.c \
	InpMtx.fullAdj.c \
	InpMtx.gmmm.c \
	InpMtx.gmvm.c \
	InpMtx.init.c \
	InpMtx.input.c \
	InpMtx.instance.c \
	InpMtx.IO.c \
	InpMtx.map.c \
	InpMtx.mvm.c \
	InpMtx.mvmVector.c \
	InpMtx.permute.c \
	InpMtx.profile.c \
	InpMtx.support.c \
	InpMtx.util.c

# headers
HDR = \
	InpMtx.h
HDR_LINK = $(HDR:.h=.h_link)

# subdirectories
SUB_DIR = 
#subdir_driver: subdir_loop
# SUB_DIR is empty
subdir_driver: 

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = $(MACRO_DIR)/empty
include $(DEPEND)
