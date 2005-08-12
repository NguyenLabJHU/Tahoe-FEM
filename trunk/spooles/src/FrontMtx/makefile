# $Id: makefile,v 1.4 2005-08-12 05:27:11 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	FrontMtx.basics.o \
	FrontMtx.factor.o \
	FrontMtx.factorUtil.o \
	FrontMtx.init.o \
	FrontMtx.instance.o \
	FrontMtx.IO.o \
	FrontMtx.loadEntries.o \
	FrontMtx.permute.o \
	FrontMtx.postponed.o \
	FrontMtx.postProcess.o \
	FrontMtx.QRfactor.o \
	FrontMtx.QRsolve.o \
	FrontMtx.QRutil.o \
	FrontMtx.solve.o \
	FrontMtx.solveUtil.o \
	FrontMtx.split.o \
	FrontMtx.storeFront.o \
	FrontMtx.update.o \
	FrontMtx.util.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = $(MACRO_DIR)/empty

# sources
SRC = \
	FrontMtx.basics.c \
	FrontMtx.factor.c \
	FrontMtx.factorUtil.c \
	FrontMtx.init.c \
	FrontMtx.instance.c \
	FrontMtx.IO.c \
	FrontMtx.loadEntries.c \
	FrontMtx.permute.c \
	FrontMtx.postponed.c \
	FrontMtx.postProcess.c \
	FrontMtx.QRfactor.c \
	FrontMtx.QRsolve.c \
	FrontMtx.QRutil.c \
	FrontMtx.solve.c \
	FrontMtx.solveUtil.c \
	FrontMtx.split.c \
	FrontMtx.storeFront.c \
	FrontMtx.update.c \
	FrontMtx.util.c

# headers
HDR = \
	FrontMtx.h
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
