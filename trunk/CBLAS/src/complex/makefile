# $Id: makefile,v 1.4 2001-02-22 23:40:23 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	caxpy.o \
	ccopy.o \
	cdotc.o \
	cgemv.o \
	cgerc.o \
	chemv.o \
	cher2.o \
	cmyblas2.o \
	cscal.o \
	ctrsv.o \
	dcabs1.o \
	icamax.o \
	scasum.o \
	scnrm2.o
OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	caxpy.c \
	ccopy.c \
	cdotc.c \
	cgemv.c \
	cgerc.c \
	chemv.c \
	cher2.c \
	cmyblas2.c \
	cscal.c \
	ctrsv.c \
	dcabs1.c \
	icamax.c \
	scasum.c \
	scnrm2.c

# headers
#HDR = 
#HDR_LINK = $(HDR:.h=.h_link)

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
