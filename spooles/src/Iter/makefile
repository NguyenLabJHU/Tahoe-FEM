# $Id: makefile,v 1.4 2005-08-12 05:27:12 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	Iter.bgmresl.o \
	Iter.bgmresr.o \
	Iter.bicgstabl.o \
	Iter.bicgstabr.o \
	Iter.DenseMtx_mmm.o \
	Iter.mlbicgstabl.o \
	Iter.mlbicgstabr.o \
	Iter.pcgl.o \
	Iter.pcgr.o \
	Iter.tfqmrl.o \
	Iter.tfqmrr.o \
	Iter.util.o \
	Iter.zbicgstabl.o \
	Iter.zbicgstabr.o \
	Iter.zmlbicgstabl.o \
	Iter.zmlbicgstabr.o \
	Iter.zpcgl.o \
	Iter.zpcgr.o \
	Iter.ztfqmrl.o \
	Iter.ztfqmrr.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = $(MACRO_DIR)/empty

# sources
SRC = \
	Iter.bgmresl.c \
	Iter.bgmresr.c \
	Iter.bicgstabl.c \
	Iter.bicgstabr.c \
	Iter.DenseMtx_mmm.c \
	Iter.mlbicgstabl.c \
	Iter.mlbicgstabr.c \
	Iter.pcgl.c \
	Iter.pcgr.c \
	Iter.tfqmrl.c \
	Iter.tfqmrr.c \
	Iter.util.c \
	Iter.zbicgstabl.c \
	Iter.zbicgstabr.c \
	Iter.zmlbicgstabl.c \
	Iter.zmlbicgstabr.c \
	Iter.zpcgl.c \
	Iter.zpcgr.c \
	Iter.ztfqmrl.c \
	Iter.ztfqmrr.c

# headers
HDR = \
	Iter.h
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
