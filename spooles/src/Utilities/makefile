# $Id: makefile,v 1.4 2005-08-12 05:27:13 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	Utilities.axpy.o \
	Utilities.CV.o \
	Utilities.DV.o \
	Utilities.FV.o \
	Utilities.iohb.o \
	Utilities.IP.o \
	Utilities.IV.o \
	Utilities.I2OP.o \
	Utilities.newsort.o \
	Utilities.PCV.o \
	Utilities.PDV.o \
	Utilities.PFV.o \
	Utilities.PIV.o \
	Utilities.sortAndCompress.o \
	Utilities.ZV.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = $(MACRO_DIR)/empty

# sources
SRC = \
	Utilities.axpy.c \
	Utilities.CV.c \
	Utilities.DV.c \
	Utilities.FV.c \
	Utilities.iohb.c \
	Utilities.IP.c \
	Utilities.IV.c \
	Utilities.I2OP.c \
	Utilities.newsort.c \
	Utilities.PCV.c \
	Utilities.PDV.c \
	Utilities.PFV.c \
	Utilities.PIV.c \
	Utilities.sortAndCompress.c \
	Utilities.ZV.c

# headers
HDR = \
	axpy_.h \
	CV_.h \
	DLIST_.h \
	DV_.h \
	DVsort_.h \
	FileType_.h \
	FV_.h \
	iohb_.h \
	IP_.h \
	IV_.h \
	IVsort_.h \
	I2OP_.h \
	mdot_.h \
	misc_.h \
	MM_.h \
	ND_.h \
	newsort_.h \
	PCV_.h \
	PDV_.h \
	PFV_.h \
	PIV_.h \
	util_.h \
	Utilities.h \
	ZV_.h
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
