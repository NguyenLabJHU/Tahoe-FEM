# $Id: makefile,v 1.4 2005-08-12 05:27:13 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	misc.drawGraphEPS.o \
	misc.fpgrids.o \
	misc.localND.o \
	misc.mkNDlinsys.o \
	misc.mkNDlinsysQR.o \
	misc.ND.o \
	misc.ND2.o \
	misc.orderViaBestOfNDandMS.o \
	misc.orderViaMMD.o \
	misc.orderViaMS.o \
	misc.orderViaND.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = $(MACRO_DIR)/empty

# sources
SRC = \
	misc.drawGraphEPS.c \
	misc.fpgrids.c \
	misc.localND.c \
	misc.mkNDlinsys.c \
	misc.mkNDlinsysQR.c \
	misc.ND.c \
	misc.ND2.c \
	misc.orderViaBestOfNDandMS.c \
	misc.orderViaMMD.c \
	misc.orderViaMS.c \
	misc.orderViaND.c

# headers
HDR = \
	misc.h
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
